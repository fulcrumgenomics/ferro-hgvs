//! Cross-test for ferro's two splice-distance ladders.
//!
//! ferro classifies intronic offsets on **two** ladders that disagree, and both
//! are correct because they answer different questions:
//!
//! | shape | type | ladder | question |
//! |---|---|---|---|
//! | A | [`SpliceSiteType`] | donor 2/6/20/50, acceptor 2/12/20/50 | ferro's fine-grained splice-distance bin |
//! | B | [`SpliceRegionBin`] | 2/8 | which Sequence Ontology term does VEP assign? |
//!
//! Shape A previously existed as **three** full copies (`reference/transcript.rs`,
//! and twice in `convert/noncoding.rs`) plus a fourth partial one
//! (`IntronicRegion::from_offset`, 2/20/50 — missing the 6- and 12-rungs), and
//! shape B as **two** (`effect/mod.rs` and the web service's CDS effect
//! handler). Each shape now has one definition and the rest derive from it.
//!
//! This module pins two things a future edit could silently break:
//!
//! 1. Every shape-A consumer still agrees with the one shape-A ladder, and
//!    every shape-B consumer with the one shape-B ladder — i.e. the
//!    de-duplication stayed de-duplicated.
//! 2. The **relationship between the two shapes** at their boundary offsets. The
//!    ladders are meant to disagree; what must not happen is one of them moving
//!    relative to the other without anyone noticing, since two surfaces of one
//!    product would then classify one variant differently for a new reason.
//!
//! The web service's `predict_cds_effect` is private and feature-gated, so it is
//! not called here; it shares `SpliceRegionBin`, which is what makes it unable
//! to drift.

use ferro_hgvs::convert::noncoding::{IntronicConsequence, IntronicRegion};
use ferro_hgvs::effect::{Consequence, EffectPredictor, SpliceRegionBin};
use ferro_hgvs::hgvs::location::{CdsPos, TxPos};
use ferro_hgvs::reference::transcript::{IntronBoundary, IntronPosition, SpliceSiteType};

/// The offsets worth pinning: every rung of both ladders, and the base either
/// side of each rung. Signed — positive is the donor side.
const BOUNDARY_OFFSETS: &[i64] = &[
    1, 2, 3, 6, 7, 8, 9, 10, 12, 13, 20, 21, 50, 51, //
    -1, -2, -3, -6, -7, -8, -9, -10, -12, -13, -20, -21, -50, -51,
];

fn cds_at(offset: i64) -> CdsPos {
    CdsPos {
        base: 100,
        offset: Some(offset),
        utr3: false,
        special: None,
    }
}

fn tx_at(offset: i64) -> TxPos {
    TxPos {
        base: 100,
        offset: Some(offset),
        downstream: false,
    }
}

fn intron_position_at(offset: i64) -> IntronPosition {
    IntronPosition {
        intron_number: 1,
        boundary: if offset > 0 {
            IntronBoundary::FivePrime
        } else {
            IntronBoundary::ThreePrime
        },
        offset,
        tx_boundary_pos: 100,
        intron_length: 1000,
    }
}

// --- Shape A: one ladder, several consumers ---------------------------------

/// The fine-grained ladder itself, rung by rung, in both directions.
#[test]
fn shape_a_bins_every_boundary_offset_as_specified() {
    use SpliceSiteType::*;
    let expected: &[(i64, SpliceSiteType)] = &[
        (1, DonorCanonical),
        (2, DonorCanonical),
        (3, DonorExtended),
        (6, DonorExtended),
        (7, DonorRegion),
        (20, DonorRegion),
        (21, NearSplice),
        (50, NearSplice),
        (51, DeepIntronic),
        (-1, AcceptorCanonical),
        (-2, AcceptorCanonical),
        (-3, AcceptorExtended),
        (-12, AcceptorExtended),
        (-13, AcceptorRegion),
        (-20, AcceptorRegion),
        (-21, NearSplice),
        (-50, NearSplice),
        (-51, DeepIntronic),
    ];
    assert!(
        !expected.is_empty(),
        "the ladder table must be non-empty, or this test passes vacuously"
    );
    for (offset, want) in expected {
        assert_eq!(
            SpliceSiteType::from_signed_offset(*offset),
            *want,
            "SpliceSiteType::from_signed_offset({offset})"
        );
    }
}

/// `IntronPosition` carries its side as a field, not as the sign of its offset,
/// so an offset of 0 on the 5' boundary must still classify as a **donor**.
///
/// This is the one case where reading the side off the sign would change
/// behaviour, which is why the shared ladder is reachable both ways
/// (`from_signed_offset` for offsets, `from_distance_on_side` for the field).
#[test]
fn shape_a_reads_its_side_from_the_boundary_field_not_the_sign() {
    let five_prime_zero = IntronPosition {
        intron_number: 1,
        boundary: IntronBoundary::FivePrime,
        offset: 0,
        tx_boundary_pos: 100,
        intron_length: 1000,
    };
    assert_eq!(
        five_prime_zero.splice_site_type(),
        SpliceSiteType::DonorCanonical,
        "a 5'-boundary position must classify as a donor even at offset 0, where \
         the sign alone would say acceptor"
    );

    // And a negatively-spelled offset on the 5' boundary is still a donor: the
    // field wins over the sign in both directions.
    let five_prime_negative = IntronPosition {
        offset: -6,
        ..five_prime_zero
    };
    assert_eq!(
        five_prime_negative.splice_site_type(),
        SpliceSiteType::DonorExtended
    );

    // The mirror: a positively-spelled offset on the 3' boundary is still an
    // acceptor, and lands on the acceptor's own 12-rung — which the donor
    // ladder does not have, so this could not pass if the side were read off
    // the sign.
    let three_prime_positive = IntronPosition {
        boundary: IntronBoundary::ThreePrime,
        offset: 12,
        ..five_prime_zero
    };
    assert_eq!(
        three_prime_positive.splice_site_type(),
        SpliceSiteType::AcceptorExtended,
        "offset +12 on the 3' boundary must be AcceptorExtended; reading the \
         side off the sign would give DonorRegion"
    );
}

/// The `+?` / `-?` unknown-offset sentinels must not panic and must not be
/// reported as a canonical splice site.
///
/// `c.<base>-?` parses to `offset == Some(i64::MIN)`
/// (`parser::position::OFFSET_UNKNOWN_NEGATIVE`), and every ladder here takes
/// the offset's magnitude. Taken with `abs()` that is a debug panic, and in
/// release it wraps back to `i64::MIN`, which satisfies every `<= n` rung and
/// reports the *closest* possible bin — a HIGH-impact canonical splice site —
/// for a position whose distance is by definition unknown. Both ladders take
/// the magnitude with `unsigned_abs`, so the sentinel lands in the farthest bin
/// instead.
///
/// This pins "does not panic, is not canonical". Whether to decline outright
/// was left open here, and #1841 settled it: the two ladders that answer *for a
/// position* — `IntronicRegion::from_offset` and the two `IntronicConsequence`
/// constructors — now return `None`, because every variant of those enums names
/// a distance band and a sentinel states that the distance is unknown. The two
/// raw bin ladders below still answer, and still answer with the farthest bin;
/// they classify a magnitude, not a position.
#[test]
fn the_unknown_offset_sentinels_do_not_panic_or_read_as_canonical() {
    for sentinel in [i64::MIN, i64::MAX] {
        assert_eq!(
            SpliceSiteType::from_signed_offset(sentinel),
            SpliceSiteType::DeepIntronic,
            "SpliceSiteType::from_signed_offset({sentinel})"
        );
        assert_eq!(
            SpliceRegionBin::from_offset(sentinel),
            SpliceRegionBin::Intron,
            "SpliceRegionBin::from_offset({sentinel})"
        );
        assert_eq!(
            IntronicRegion::from_offset(sentinel),
            None,
            "IntronicRegion::from_offset({sentinel})"
        );
        assert_eq!(
            IntronicConsequence::from_cds_pos(&cds_at(sentinel)),
            None,
            "IntronicConsequence::from_cds_pos({sentinel})"
        );
        assert_eq!(
            IntronicConsequence::from_tx_pos(&tx_at(sentinel)),
            None,
            "IntronicConsequence::from_tx_pos({sentinel})"
        );
        assert!(
            EffectPredictor::new()
                .classify_splice_variant(sentinel)
                .is_none(),
            "EffectPredictor::classify_splice_variant({sentinel})"
        );

        // `IntronPosition`'s four standalone predicates take the same magnitude
        // and had the same `abs()` hazard.
        let pos = intron_position_at(sentinel);
        assert!(pos.is_deep_intronic());
        assert!(!pos.is_canonical_splice_site());
        assert!(!pos.is_near_splice_site());
        assert!(!pos.is_extended_splice_region());
    }
}

/// `IntronicConsequence` has three entry points; all three must land on the one
/// shape-A ladder rather than on a private copy of its thresholds.
#[test]
fn every_shape_a_consumer_agrees_with_the_one_shape_a_ladder() {
    let mut checked = 0usize;
    for &offset in BOUNDARY_OFFSETS {
        let site = SpliceSiteType::from_signed_offset(offset);
        let want = IntronicConsequence::from_splice_site_type(site);

        assert_eq!(
            IntronicConsequence::from_cds_pos(&cds_at(offset)),
            Some(want),
            "IntronicConsequence::from_cds_pos at offset {offset}"
        );
        assert_eq!(
            IntronicConsequence::from_tx_pos(&tx_at(offset)),
            Some(want),
            "IntronicConsequence::from_tx_pos at offset {offset}"
        );
        assert_eq!(
            IntronicConsequence::from_intron_position(&intron_position_at(offset)),
            want,
            "IntronicConsequence::from_intron_position at offset {offset}"
        );
        checked += 1;
    }
    assert!(
        checked > 0,
        "VACUOUS: no offset was checked against the shape-A ladder"
    );
}

/// The projection `IntronicConsequence::from_splice_site_type` performs, pinned
/// against literals rather than against itself.
///
/// The test above derives its expectation *from* this projection, so it pins
/// that the three entry points agree with each other and says nothing about
/// whether the projection is right. Since de-duplicating those three onto it is
/// this change's only functional edit, the mapping needs its own outright pin —
/// the same treatment `IntronicRegion`'s historical 2/20/50 spelling gets below.
#[test]
fn the_consequence_projection_maps_every_bin_as_specified() {
    use IntronicConsequence::*;
    use SpliceSiteType::*;
    let expected: &[(SpliceSiteType, IntronicConsequence)] = &[
        (DonorCanonical, SpliceDonorVariant),
        (AcceptorCanonical, SpliceAcceptorVariant),
        (DonorExtended, SpliceDonorRegionVariant),
        (AcceptorExtended, SpliceAcceptorRegionVariant),
        (DonorRegion, SpliceRegionVariant),
        (AcceptorRegion, SpliceRegionVariant),
        (NearSplice, NearSpliceSiteVariant),
        (DeepIntronic, IntronVariant),
    ];
    assert!(
        !expected.is_empty(),
        "the projection table must be non-empty, or this test passes vacuously"
    );
    for (site, want) in expected {
        assert_eq!(
            IntronicConsequence::from_splice_site_type(*site),
            *want,
            "IntronicConsequence::from_splice_site_type({site:?})"
        );
    }
}

/// `IntronicRegion` is the side-agnostic collapse of shape A. It used to restate
/// a *partial* copy of the ladder (2/20/50), so this pins that the rungs it
/// drops are dropped by the mapping and not by a second editable ladder.
#[test]
fn the_coarse_region_view_is_a_collapse_of_the_same_ladder() {
    use IntronicRegion::*;
    let mut checked = 0usize;
    for &offset in BOUNDARY_OFFSETS {
        let want = match SpliceSiteType::from_signed_offset(offset) {
            SpliceSiteType::DonorCanonical | SpliceSiteType::AcceptorCanonical => {
                CanonicalSpliceSite
            }
            SpliceSiteType::DonorExtended
            | SpliceSiteType::AcceptorExtended
            | SpliceSiteType::DonorRegion
            | SpliceSiteType::AcceptorRegion => ExtendedSpliceRegion,
            SpliceSiteType::NearSplice => NearSpliceSite,
            SpliceSiteType::DeepIntronic => DeepIntronic,
        };
        assert_eq!(
            IntronicRegion::from_offset(offset),
            Some(want),
            "IntronicRegion::from_offset({offset})"
        );
        checked += 1;
    }
    assert!(checked > 0, "VACUOUS: no offset was checked");

    // The historical 2/20/50 spelling, pinned outright so the collapse cannot
    // quietly change what this view reports.
    assert_eq!(IntronicRegion::from_offset(2), Some(CanonicalSpliceSite));
    assert_eq!(IntronicRegion::from_offset(3), Some(ExtendedSpliceRegion));
    assert_eq!(IntronicRegion::from_offset(20), Some(ExtendedSpliceRegion));
    assert_eq!(IntronicRegion::from_offset(21), Some(NearSpliceSite));
    assert_eq!(IntronicRegion::from_offset(50), Some(NearSpliceSite));
    assert_eq!(IntronicRegion::from_offset(51), Some(DeepIntronic));
}

// --- Shape B: one ladder, several consumers ---------------------------------

#[test]
fn shape_b_bins_every_boundary_offset_as_specified() {
    use SpliceRegionBin::*;
    let expected: &[(i64, SpliceRegionBin)] = &[
        (1, Canonical),
        (2, Canonical),
        (3, Region),
        (8, Region),
        (9, Intron),
        (51, Intron),
        (-2, Canonical),
        (-8, Region),
        (-9, Intron),
    ];
    assert!(
        !expected.is_empty(),
        "the ladder table must be non-empty, or this test passes vacuously"
    );
    for (offset, want) in expected {
        assert_eq!(
            SpliceRegionBin::from_offset(*offset),
            *want,
            "SpliceRegionBin::from_offset({offset})"
        );
    }
}

#[test]
fn the_effect_predictor_agrees_with_the_one_shape_b_ladder() {
    let predictor = EffectPredictor::new();
    let mut checked = 0usize;
    for &offset in BOUNDARY_OFFSETS {
        let want = match SpliceRegionBin::from_offset(offset) {
            SpliceRegionBin::Canonical => {
                if offset > 0 {
                    Consequence::SpliceDonorVariant
                } else {
                    Consequence::SpliceAcceptorVariant
                }
            }
            SpliceRegionBin::Region => Consequence::SpliceRegionVariant,
            SpliceRegionBin::Intron => Consequence::IntronVariant,
        };
        let effect = predictor
            .classify_splice_variant(offset)
            .expect("a measured offset is not a sentinel, so it classifies");
        assert_eq!(
            effect.consequences,
            vec![want],
            "EffectPredictor::classify_splice_variant({offset})"
        );
        assert_eq!(effect.intronic_offset, Some(offset));
        checked += 1;
    }
    assert!(
        checked > 0,
        "VACUOUS: no offset was checked against the shape-B ladder"
    );
}

// --- The relationship between the two shapes --------------------------------

/// The two ladders disagree **by design**, and this pins a representative
/// sample of where.
///
/// Deliberately a sample and not an exhaustive set — the two disagree over
/// roughly `3..=50` in both directions, and pinning every one of those offsets
/// would restate a ladder rather than pin a relationship. What is pinned is one
/// offset per *kind* of disagreement, which is what an edit to either ladder
/// would move.
///
/// If a future edit makes this test fail, the question to answer is not "which
/// ladder is wrong" — it is "did I mean to change what one of ferro's two
/// surfaces reports for this offset?". Neither shape is a bug.
#[test]
fn the_two_shapes_disagree_at_the_documented_offsets() {
    // (offset, shape A, shape B) at representative offsets where the two
    // ladders put the position in bins of different coarseness. Offsets are
    // signed: positive is the donor side.
    let disagreements: &[(i64, SpliceSiteType, SpliceRegionBin)] = &[
        // +3..=+6 donor: A says "extended donor consensus", B says "splice region".
        (6, SpliceSiteType::DonorExtended, SpliceRegionBin::Region),
        // Beyond +8/-8, B has already dropped to a plain intron variant while A
        // is still inside a splice bin out to 20. This is the widest gap between
        // the two. Note A's bin differs by side here: +10 is `DonorRegion`
        // (donor's 6-rung is behind it), while -10 is still `AcceptorExtended`
        // (the acceptor's extended rung reaches 12).
        (10, SpliceSiteType::DonorRegion, SpliceRegionBin::Intron),
        (
            -10,
            SpliceSiteType::AcceptorExtended,
            SpliceRegionBin::Intron,
        ),
        (
            -12,
            SpliceSiteType::AcceptorExtended,
            SpliceRegionBin::Intron,
        ),
        (-20, SpliceSiteType::AcceptorRegion, SpliceRegionBin::Intron),
    ];
    assert!(
        !disagreements.is_empty(),
        "VACUOUS: the two shapes must be pinned as disagreeing somewhere, or \
         this guard would pass if they were silently merged"
    );
    for (offset, want_a, want_b) in disagreements {
        assert_eq!(
            SpliceSiteType::from_signed_offset(*offset),
            *want_a,
            "shape A at offset {offset}"
        );
        assert_eq!(
            SpliceRegionBin::from_offset(*offset),
            *want_b,
            "shape B at offset {offset}"
        );
    }
}

/// The converse half: the two shapes must still agree at the canonical site and
/// deep in the intron. A change that made them disagree *there* would be a real
/// defect rather than the designed difference above.
#[test]
fn the_two_shapes_agree_at_the_canonical_site_and_deep_intronic() {
    for offset in [1i64, 2, -1, -2] {
        assert!(
            matches!(
                SpliceSiteType::from_signed_offset(offset),
                SpliceSiteType::DonorCanonical | SpliceSiteType::AcceptorCanonical
            ),
            "shape A must call offset {offset} canonical"
        );
        assert_eq!(
            SpliceRegionBin::from_offset(offset),
            SpliceRegionBin::Canonical
        );
    }
    for offset in [51i64, 200, -51, -200] {
        assert_eq!(
            SpliceSiteType::from_signed_offset(offset),
            SpliceSiteType::DeepIntronic,
            "shape A must call offset {offset} deep intronic"
        );
        assert_eq!(
            SpliceRegionBin::from_offset(offset),
            SpliceRegionBin::Intron
        );
    }
}
