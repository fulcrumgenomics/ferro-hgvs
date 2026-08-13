//! Regression guards for #1767 — the `+?` / `-?` unknown-offset sentinels must
//! never be classified as if they were a *measured* intronic distance.
//!
//! `parse_offset` represents `+?` as [`OFFSET_UNKNOWN_POSITIVE`] (`i64::MAX`)
//! and `-?` as [`OFFSET_UNKNOWN_NEGATIVE`] (`i64::MIN`), and `CdsPos::is_intronic`
//! is satisfied by both. So a description spelled `…:c.100-?del` carried the
//! sentinel straight into every splice classifier's `.abs()`:
//!
//! * **debug** — `i64::MIN.abs()` is `attempt to negate with overflow`, a panic;
//! * **release** — `i64::MIN.wrapping_abs()` is `i64::MIN`, which is *negative*,
//!   so every `abs_offset <= 2` test passed and an unknown offset was reported
//!   as a canonical splice site with HIGH impact.
//!
//! The two profiles disagreed, which is why the assertions below are written
//! against the *value* rather than against a panic: after the fix both profiles
//! produce the same answer, so one set of assertions covers both. The debug
//! panic is covered implicitly — these tests run under `overflow-checks`.
//!
//! The contract applied here is the one `CdsPos::has_unknown_offset` already
//! documents (closing #1087): the sentinels denote an unknown position
//! unbounded in one direction, so no distance can be derived from them at all.
//! **Every** classifier now says so by returning `None` — #1767 could only give
//! that answer where the return type already admitted it, and #1841 took the
//! public-API break for the two that did not (`IntronicRegion::from_offset` and
//! `EffectPredictor::classify_splice_variant`). The guards that used to pin the
//! fall-back "intronic" answer for those two now pin the decline; see
//! `tests/it/issue_1841_option_returning_splice_classifiers.rs` for the
//! break's own guards.
//!
//! [`OFFSET_UNKNOWN_POSITIVE`]: ferro_hgvs::hgvs::parser::position::OFFSET_UNKNOWN_POSITIVE
//! [`OFFSET_UNKNOWN_NEGATIVE`]: ferro_hgvs::hgvs::parser::position::OFFSET_UNKNOWN_NEGATIVE

use ferro_hgvs::convert::coding::validate_cds_pos;
use ferro_hgvs::convert::noncoding::{
    is_clinically_significant_splice_position, validate_tx_pos, IntronicConsequence, IntronicRegion,
};
use ferro_hgvs::effect::{Consequence, EffectPredictor, Impact};
use ferro_hgvs::hgvs::location::{CdsPos, IvsNotation, TxPos};
use ferro_hgvs::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::transcript::{
    Exon, GenomeBuild, IntronBoundary, IntronPosition, ManeStatus, SpliceSiteType, Strand,
    Transcript,
};

/// Both unknown-offset sentinels, so every guard below is exercised
/// symmetrically. #1767 notes that a fix teaching the classifiers only about
/// `i64::MIN` would leave `+?` classifying as deep-intronic "by accident
/// rather than by decision".
const SENTINELS: [i64; 2] = [OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE];

fn cds_with_offset(offset: i64) -> CdsPos {
    CdsPos {
        base: 100,
        offset: Some(offset),
        utr3: false,
        special: None,
    }
}

fn tx_with_offset(offset: i64) -> TxPos {
    TxPos {
        base: 100,
        offset: Some(offset),
        downstream: false,
    }
}

/// A minimal coding transcript, built programmatically, whose only job is to
/// let the two validators get as far as the offset check. Base 100 is well
/// inside its CDS, so anything the validators reject is about the offset.
fn single_exon_transcript() -> Transcript {
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        "A".repeat(300),
        Some(1),
        Some(300),
        vec![Exon::new(1, 1, 300)],
        None,
        None,
        None,
        GenomeBuild::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

/// The premise of the issue: `c.100-?del` really does parse to the sentinel,
/// and really does satisfy `is_intronic`, so it reaches the classifiers.
/// Asserted against the AST rather than the rendered string — a round-trip
/// would hold even if the offset were carried some other way.
#[test]
fn an_unknown_offset_description_carries_the_sentinel_into_the_classifiers() {
    let variant = parse_hgvs("NM_000088.3:c.100-?del").expect("`c.100-?del` is legal HGVS");
    let HgvsVariant::Cds(cds) = &variant else {
        panic!("expected a c. variant, got {variant:?}");
    };
    let start = cds
        .loc_edit
        .location
        .start
        .inner()
        .expect("`c.100-?` has a certain start position");

    assert_eq!(
        start.offset,
        Some(OFFSET_UNKNOWN_NEGATIVE),
        "`-?` must parse to the sentinel — the whole issue rests on this"
    );
    assert!(
        start.is_intronic(),
        "`is_intronic` is `offset.is_some() && offset != Some(0)`, so the \
         sentinel satisfies it and reaches every splice classifier"
    );
    assert!(
        start.has_unknown_offset(),
        "and the shared predicate recognises it, which is what the fix keys on"
    );
    assert!(
        variant.to_string().contains("100-?"),
        "the sentinel still renders as `-?`, so no output moves"
    );
}

#[test]
fn intronic_consequence_declines_an_unknown_cds_offset() {
    for offset in SENTINELS {
        assert_eq!(
            IntronicConsequence::from_cds_pos(&cds_with_offset(offset)),
            None,
            "an unknown offset ({offset}) has no derivable distance, so it cannot \
             be classified as a splice consequence"
        );
    }
}

#[test]
fn intronic_consequence_declines_an_unknown_tx_offset() {
    for offset in SENTINELS {
        assert_eq!(
            IntronicConsequence::from_tx_pos(&tx_with_offset(offset)),
            None,
            "an unknown offset ({offset}) has no derivable distance on the n. axis either"
        );
    }
}

#[test]
fn intronic_region_declines_an_unknown_offset() {
    for offset in SENTINELS {
        assert_eq!(
            IntronicRegion::from_cds_pos(&cds_with_offset(offset)),
            None,
            "every `IntronicRegion` variant asserts a distance band, so an unknown \
             offset ({offset}) must decline rather than pick one"
        );
        assert_eq!(
            IntronicRegion::from_tx_pos(&tx_with_offset(offset)),
            None,
            "same on the n. axis for offset {offset}"
        );
    }
}

/// The release-mode outcome the issue calls out: an unknown offset wrapped to a
/// negative `abs`, passed `<= 2`, and was reported as clinically significant.
#[test]
fn an_unknown_offset_is_not_a_clinically_significant_splice_position() {
    for offset in SENTINELS {
        assert!(
            !is_clinically_significant_splice_position(&cds_with_offset(offset)),
            "an unknown offset ({offset}) is not *known* to affect splicing"
        );
    }
}

/// The property #1767 was filed for, stated in the form it takes after #1841
/// widened the answer from "the weakest true statement" to "no statement".
///
/// What must never happen is a *splice-site* claim — that was the release-mode
/// defect (`SpliceAcceptorVariant`/HIGH). A decline satisfies that strictly more
/// than `IntronVariant` did, so this guard is not weakened by the change; it is
/// the same claim asserted against a return type that can carry it.
#[test]
fn classify_splice_variant_makes_no_distance_claim_for_an_unknown_offset() {
    let predictor = EffectPredictor::new();
    for offset in SENTINELS {
        assert!(
            predictor.classify_splice_variant(offset).is_none(),
            "an unknown offset ({offset}) must not be reported as a splice \
             donor/acceptor or splice-region variant, and since #1841 it is not \
             reported as anything"
        );
    }
}

/// The #1765 lesson applied here: a guard keyed on sentinel *identity* leaves
/// the neighbouring magnitudes broken. `i64::MIN` is the only value whose
/// `.abs()` overflows, so the classifiers must be total for every `i64` —
/// not merely for the two values that happen to be spelled `?`.
#[test]
fn the_classifiers_are_total_over_every_offset_magnitude() {
    let predictor = EffectPredictor::new();
    for offset in [i64::MIN, i64::MIN + 1, i64::MAX - 1, i64::MAX, 0, -1, 1] {
        // Each of these overflowed (debug: panic) or wrapped (release) before
        // the fix; the assertion is simply that they now return.
        let _ = predictor.classify_splice_variant(offset);
        let _ = IntronicConsequence::from_cds_pos(&cds_with_offset(offset));
        let _ = IntronicConsequence::from_tx_pos(&tx_with_offset(offset));
        let _ = IntronicRegion::from_cds_pos(&cds_with_offset(offset));
        let _ = IntronicRegion::from_tx_pos(&tx_with_offset(offset));
        let _ = is_clinically_significant_splice_position(&cds_with_offset(offset));
        let _ = IntronicRegion::from_offset(offset);
    }
}

/// `IntronPosition` is normally derived from genomic coordinates, but its
/// fields are public, and its five predicates all read `offset` as a magnitude.
/// They must be total: before the fix each overflowed on `i64::MIN`, and in
/// release the wrapped negative value made `is_canonical_splice_site` answer
/// `true` for an offset that states no distance.
#[test]
fn the_intron_position_predicates_are_total() {
    for offset in [i64::MIN, i64::MIN + 1, i64::MAX - 1, i64::MAX] {
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::ThreePrime,
            offset,
            tx_boundary_pos: 100,
            intron_length: 1_000,
        };
        assert!(
            pos.is_deep_intronic(),
            "a magnitude of {offset} is not within 50bp of an exon"
        );
        assert!(!pos.is_canonical_splice_site(), "offset {offset}");
        assert!(!pos.is_near_splice_site(), "offset {offset}");
        assert!(!pos.is_extended_splice_region(), "offset {offset}");
        assert_eq!(
            pos.splice_site_type(),
            SpliceSiteType::DeepIntronic,
            "offset {offset}"
        );
    }
}

/// `IvsPos` carries the same two predicate names and the same defect, and
/// `IvsNotation::to_ivs` maps a `CdsPos`/`TxPos` offset straight into it — so the
/// sentinel reaches it from exactly the description this issue is about.
#[test]
fn the_ivs_predicates_are_total_and_reachable_from_an_unknown_offset() {
    for offset in SENTINELS {
        let ivs = cds_with_offset(offset)
            .to_ivs(1)
            .expect("a position with an offset converts to IVS notation");
        assert_eq!(ivs.offset, offset, "the offset is carried through verbatim");
        assert!(!ivs.is_canonical_splice_site(), "offset {offset}");
        assert!(ivs.is_deep_intronic(), "offset {offset}");
    }
}

/// Coordinate conversion must decline an unknown offset distinguishably rather
/// than overflow on it. Both validators previously panicked under
/// `overflow-checks` for `-?`; `+?` was rejected, but as an out-of-range
/// number rather than as the unknown it is.
#[test]
fn the_coordinate_validators_decline_an_unknown_offset() {
    let transcript = single_exon_transcript();

    for offset in SENTINELS {
        let cds_err = validate_cds_pos(&cds_with_offset(offset), &transcript)
            .expect_err("an unknown offset denotes no coordinate");
        assert!(
            cds_err.to_string().contains("Unknown intronic offset"),
            "offset {offset} must be declined as unknown, not as out of range; got {cds_err}"
        );

        let tx_err = validate_tx_pos(&tx_with_offset(offset), &transcript)
            .expect_err("an unknown offset denotes no coordinate");
        assert!(
            tx_err.to_string().contains("Unknown intronic offset"),
            "offset {offset} on the n. axis; got {tx_err}"
        );
    }

    // A measured intronic offset still validates, so the guard is scoped to
    // the sentinels rather than rejecting intronic positions generally.
    assert!(validate_cds_pos(&cds_with_offset(-2), &transcript).is_ok());
    assert!(validate_tx_pos(&tx_with_offset(-2), &transcript).is_ok());
}

/// Guard the fix against over-reach: *measured* offsets must classify exactly
/// as they did before. These are the values the existing unit tests pin.
#[test]
fn measured_offsets_are_unaffected() {
    let predictor = EffectPredictor::new();

    // `.expect` rather than a bare unwrap: since #1841 a `None` here would mean
    // the decline had widened past the sentinels, which is the failure worth
    // naming.
    let classify = |offset: i64| {
        predictor
            .classify_splice_variant(offset)
            .unwrap_or_else(|| panic!("a measured offset ({offset}) must still classify"))
    };

    // Impact is asserted alongside the consequence because it is the half the
    // #1767 defect got wrong (HIGH for a position stating no distance), and
    // because this file's only previous `Impact` assertion was on the sentinel
    // — which since #1841 has no effect to carry one.
    for (offset, consequence, impact) in [
        (1, Consequence::SpliceDonorVariant, Impact::High),
        (-2, Consequence::SpliceAcceptorVariant, Impact::High),
        (5, Consequence::SpliceRegionVariant, Impact::Low),
        (50, Consequence::IntronVariant, Impact::Modifier),
    ] {
        let effect = classify(offset);
        assert_eq!(effect.consequences, vec![consequence]);
        assert_eq!(effect.impact, impact);
    }

    assert_eq!(
        IntronicConsequence::from_cds_pos(&cds_with_offset(-2)),
        Some(IntronicConsequence::SpliceAcceptorVariant)
    );
    assert_eq!(
        IntronicConsequence::from_cds_pos(&cds_with_offset(1)),
        Some(IntronicConsequence::SpliceDonorVariant)
    );
    assert_eq!(
        IntronicRegion::from_cds_pos(&cds_with_offset(-2)),
        Some(IntronicRegion::CanonicalSpliceSite)
    );
    assert_eq!(
        IntronicRegion::from_cds_pos(&cds_with_offset(100)),
        Some(IntronicRegion::DeepIntronic)
    );
    assert!(is_clinically_significant_splice_position(&cds_with_offset(
        -2
    )));
}
