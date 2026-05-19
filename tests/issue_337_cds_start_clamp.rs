//! Regression test for issue #337: 5'-direction shuffle must clamp at
//! CDS-start so a CDS variant does not silently shift into the 5'UTR
//! (changing both the position and the coordinate sub-axis from
//! `c.<N>` to `c.-N`).
//!
//! Biocommons case (from #324 baseline-failures):
//!   `NM_001001656.1:c.1del`  (5prime+cross, 5prime+no-cross)
//!     expected: c.1del
//!     ferro:    c.-2del   (shifts through a 5'UTR A-tract)
//!
//! The synthetic transcript here uses 5'UTR `AAA` immediately before a
//! CDS that also starts with `A`, so a 5'-direction shuffle of `c.1del`
//! can reach `c.-3del` if unclamped — matching the biocommons failure
//! pattern.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript:
///   tx 1-3   : 5'UTR `AAA`   (c.-3 .. c.-1)
///   tx 4-12  : CDS   `ATGAAATAG`  (c.1 .. c.9 — start codon, AAA tract, stop)
///   tx 13-20 : 3'UTR `CCCCCCCC`   (c.*1 .. c.*8)
///
/// CDS starts at tx 4 (`A` of `ATG`), so `c.1` = `A`. Walking 5' from
/// `c.1`, the immediate upstream bases are `A A A` (5'UTR) — without
/// the CDS-start clamp, a 5'-direction shuffle of `c.1del` shifts to
/// `c.-3del`.
fn provider_with_utr_padded_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "AAAATGAAATAGCCCCCCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(4),
        Some(12),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

fn normalize_with(direction: ShuffleDirection, cross: bool, input: &str) -> String {
    let mut config = NormalizeConfig::default().with_direction(direction);
    if cross {
        config = config.allow_crossing_boundaries();
    }
    let normalizer = Normalizer::with_config(provider_with_utr_padded_transcript(), config);
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn five_prime_no_cross_does_not_shift_c1_del_into_5utr() {
    // c.1del under 5prime + cross_boundaries=false must stay at c.1del
    // — the 5'-direction shuffle clamps at the CDS-start regardless of
    // upstream homopolymers in the 5'UTR.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, false, "NM_TEST.1:c.1del"),
        "NM_TEST.1:c.1del",
    );
}

#[test]
fn five_prime_cross_does_not_shift_c1_del_into_5utr() {
    // Same input under cross_boundaries=true. cross_boundaries authorizes
    // crossing exon-intron boundaries, NOT the CDS↔UTR sub-axis boundary.
    // Per #337's spec note, 5'-direction shifting across c.1 lands in c.-N
    // which is a different coordinate sub-axis; the clamp applies in
    // both cross modes.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, true, "NM_TEST.1:c.1del"),
        "NM_TEST.1:c.1del",
    );
}

#[test]
fn five_prime_no_cross_does_not_shift_c1_dup_into_5utr() {
    // c.1dup under 5prime + cross_boundaries=false: same shape as the
    // c.1del case for NM_001001656.1:c.1dup in the biocommons corpus.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, false, "NM_TEST.1:c.1dup"),
        "NM_TEST.1:c.1dup",
    );
}

#[test]
fn five_prime_cross_does_not_shift_c1_dup_into_5utr() {
    // c.1dup × cross=true — completes the (del × dup) × (cross={true,false})
    // matrix called out by the biocommons corpus and the review.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, true, "NM_TEST.1:c.1dup"),
        "NM_TEST.1:c.1dup",
    );
}

#[test]
fn three_prime_shift_inside_cds_still_shuffles() {
    // Sanity check that we didn't break the 3'-direction path. The CDS
    // is `ATGAAATAG`; c.4del deletes one `A` of the c.4_c.6 `AAA` tract
    // so 3'-direction shift should move it to c.6del.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, false, "NM_TEST.1:c.4del"),
        "NM_TEST.1:c.6del",
    );
}

#[test]
fn five_prime_inside_cds_still_shuffles() {
    // c.6del inside the CDS-internal AAA tract; under 5prime+no-cross,
    // the deletion should shift 5' but only within CDS — landing at c.4.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, false, "NM_TEST.1:c.6del"),
        "NM_TEST.1:c.4del",
    );
}

/// Follow-up #350: a variant whose start and end map to different
/// coordinate sub-axes (5'UTR / CDS / 3'UTR) has no well-defined 3'-rule
/// shuffle, so `normalize_cds` short-circuits — preserving the input
/// position verbatim — and emits `NormalizationWarning::CrossAxisVariantNotShuffled`.
#[test]
fn cross_axis_5utr_to_cds_does_not_shuffle_and_warns() {
    // `c.-1_1del` against cds_start=4: tx_start=3 (5'UTR), tx_end=4 (CDS).
    // The straddling range crosses the 5'UTR↔CDS axis, so the result is
    // returned unchanged (modulo canonicalization) with a warning.
    //
    // Asserted under both cross_boundaries modes: the cross-axis bail is a
    // sub-axis (CDS↔UTR) check, orthogonal to exon-crossing policy, so the
    // warning must fire regardless of `allow_crossing_boundaries`.
    for cross in [false, true] {
        let mut config = NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime);
        if cross {
            config = config.allow_crossing_boundaries();
        }
        let normalizer = Normalizer::with_config(provider_with_utr_padded_transcript(), config);
        let variant = parse_hgvs("NM_TEST.1:c.-1_1del").expect("parse");
        let result = normalizer
            .normalize_with_warnings(&variant)
            .expect("normalize");
        assert_eq!(
            format!("{}", result.result),
            "NM_TEST.1:c.-1_1del",
            "cross={cross}",
        );
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.code() == "CROSS_AXIS_VARIANT_NOT_SHUFFLED"),
            "cross={cross}: expected CROSS_AXIS_VARIANT_NOT_SHUFFLED warning, got {:?}",
            result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
        );
    }
}

#[test]
fn cross_axis_cds_to_3utr_does_not_shuffle_and_warns() {
    // `c.9_*1del` against cds_end=12: tx_start=12 (CDS), tx_end=13 (3'UTR).
    // Same bail behavior across the CDS↔3'UTR boundary, under both
    // cross_boundaries modes.
    for cross in [false, true] {
        let mut config = NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime);
        if cross {
            config = config.allow_crossing_boundaries();
        }
        let normalizer = Normalizer::with_config(provider_with_utr_padded_transcript(), config);
        let variant = parse_hgvs("NM_TEST.1:c.9_*1del").expect("parse");
        let result = normalizer
            .normalize_with_warnings(&variant)
            .expect("normalize");
        assert_eq!(
            format!("{}", result.result),
            "NM_TEST.1:c.9_*1del",
            "cross={cross}",
        );
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.code() == "CROSS_AXIS_VARIANT_NOT_SHUFFLED"),
            "cross={cross}: expected CROSS_AXIS_VARIANT_NOT_SHUFFLED warning, got {:?}",
            result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
        );
    }
}

#[test]
fn within_axis_cds_variant_does_not_emit_cross_axis_warning() {
    // Sanity check: a purely-CDS variant must NOT emit the warning.
    let normalizer = Normalizer::with_config(
        provider_with_utr_padded_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_TEST.1:c.1del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert!(
        result
            .warnings
            .iter()
            .all(|w| w.code() != "CROSS_AXIS_VARIANT_NOT_SHUFFLED"),
        "unexpected CROSS_AXIS_VARIANT_NOT_SHUFFLED warning on within-axis variant",
    );
}

/// Follow-up #349: when the axis clamp constrains the 3'-rule shuffle
/// tighter than the exon bound would allow, `normalize_cds` emits
/// `NormalizationWarning::AxisClampApplied` so callers can see why the
/// canonical position did not shift further.
#[test]
fn axis_clamp_emits_warning_when_constraining_5prime_shuffle() {
    // c.1del under 5prime+no-cross. Without the axis clamp it would shift
    // back into the 5'UTR (per fixture: 5'UTR is `AAA`, CDS starts with
    // `A`). The clamp holds the result at c.1del. The clamp firing must
    // emit AxisClampApplied so callers can flag for review.
    let normalizer = Normalizer::with_config(
        provider_with_utr_padded_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_TEST.1:c.1del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_TEST.1:c.1del");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "AXIS_CLAMP_APPLIED"),
        "expected AXIS_CLAMP_APPLIED warning, got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
    );
}

#[test]
fn no_axis_clamp_when_upstream_utr_base_does_not_match() {
    // Pins the false-positive class CodeRabbit flagged for #343:
    // `NM_INTRA.1:c.1del` under 5prime + cross_boundaries=true on a
    // fixture whose upstream UTR base does NOT match the deleted CDS
    // base. The clamped shuffle and the unclamped exon-only shuffle
    // both stop at c.1del because the upstream UTR `C` doesn't match
    // the deleted `A`, so the axis clamp never actually constrained
    // movement. The warning must NOT fire.
    //
    // NM_INTRA.1 layout (defined in the homopolymer test below):
    //   5'UTR `CCC` (non-A), CDS `AAAATGTAG`, 3'UTR `CCCCCCCC`.
    let mut provider = MockProvider::new();
    let sequence = "CCCAAAATGTAGCCCCCCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_INTRA.1".to_string(),
        Some("INTRA".to_string()),
        Strand::Plus,
        sequence,
        Some(4),
        Some(12),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default()
            .with_direction(ShuffleDirection::FivePrime)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs("NM_INTRA.1:c.1del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_INTRA.1:c.1del");
    assert!(
        result
            .warnings
            .iter()
            .all(|w| w.code() != "AXIS_CLAMP_APPLIED"),
        "unexpected AXIS_CLAMP_APPLIED — upstream UTR base ('C') doesn't \
         match deleted CDS base ('A'), so the shuffle never tried to \
         cross the boundary and the clamp was inert. Warnings: {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
}

#[test]
fn no_axis_clamp_when_5prime_shuffle_stays_within_axis() {
    // c.6del 5prime shifts to c.4del — entirely within CDS, no clamp.
    let normalizer = Normalizer::with_config(
        provider_with_utr_padded_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_TEST.1:c.6del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_TEST.1:c.4del");
    assert!(
        result
            .warnings
            .iter()
            .all(|w| w.code() != "AXIS_CLAMP_APPLIED"),
        "unexpected AXIS_CLAMP_APPLIED on shuffle that stays within axis",
    );
}

#[test]
fn five_prime_shuffles_intra_cds_homopolymer_all_the_way_to_c1() {
    // Probe for the off-by-one Opus flagged in review. Build a fresh
    // transcript whose CDS *starts* with `AAA...` so a 5'-direction
    // shuffle of `c.3del` would (correctly) walk back to `c.1del`.
    // If the axis clamp is one HGVS position too restrictive, ferro
    // would stop at `c.2del` instead.
    //
    // Transcript layout:
    //   tx 1-3   : 5'UTR `CCC`  (non-A, so the upstream-of-CDS predicate fails)
    //   tx 4-12  : CDS   `AAAATGTAG` (starts with `AAAA`, c.1..c.4 all = `A`)
    //   tx 13-20 : 3'UTR `CCCCCCCC`
    let mut provider = MockProvider::new();
    let sequence = "CCCAAAATGTAGCCCCCCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_INTRA.1".to_string(),
        Some("INTRA".to_string()),
        Strand::Plus,
        sequence,
        Some(4),
        Some(12),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);

    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_INTRA.1:c.3del").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(
        format!("{}", normalized),
        "NM_INTRA.1:c.1del",
        "5'-shuffle through an intra-CDS A-tract must reach c.1del, not stop at c.2del",
    );
}
