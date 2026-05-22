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

/// Intra-CDS-homopolymer fixture: a CDS whose first 4 bases are `AAAA`
/// (so the 5'UTR/CDS axis clamp must let a 5'-direction `c.3del` walk
/// inside the CDS A-tract all the way to `c.1del`), and whose 5'UTR
/// uses non-A bases so the upstream-of-CDS predicate fails. Used by
/// `five_prime_shuffles_intra_cds_homopolymer_all_the_way_to_c1` and
/// `no_axis_clamp_when_upstream_utr_base_does_not_match`.
///
/// Transcript layout:
///   tx 1-3   : 5'UTR `CCC`  (non-A — the clamp is inert because the
///                           upstream UTR base does not match the
///                           deleted CDS base)
///   tx 4-12  : CDS   `AAAATGTAG` (starts with `AAAA`, c.1..c.4 all `A`)
///   tx 13-20 : 3'UTR `CCCCCCCC`
fn provider_with_intra_cds_homopolymer_transcript() -> MockProvider {
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
    let normalizer = Normalizer::with_config(
        provider_with_intra_cds_homopolymer_transcript(),
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
    // Probe for the off-by-one Opus flagged in review. The fixture
    // `provider_with_intra_cds_homopolymer_transcript`'s CDS *starts*
    // with `AAA...` so a 5'-direction shuffle of `c.3del` would
    // (correctly) walk back to `c.1del`. If the axis clamp is one HGVS
    // position too restrictive, ferro would stop at `c.2del` instead.
    let normalizer = Normalizer::with_config(
        provider_with_intra_cds_homopolymer_transcript(),
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

// =============================================================================
// PR #343 deferred follow-ups (issue #391 item 2)
//
// PR #343 (CDS↔UTR axis clamp for the 3'-rule shuffle) deferred two
// symmetry tests in its body:
//
//   - "5'-direction `c.*N` clamp at `c.*1` symmetry test"
//   - "minus-strand sanity test"
//
// Both pin behavior that ALREADY works on the existing clamp logic —
// they were called out as gaps in the test matrix, not as missing
// fixes. This block lands them, plus the natural 3'-direction CDS-end
// symmetry (a `c.<N>` variant at `c.cds_end` trying to 3'-shift INTO
// the 3'UTR), so a future axis-clamp refactor cannot regress any
// quadrant of the (direction × axis × strand × edit-kind) matrix.
// =============================================================================

/// Mirror fixture of `provider_with_utr_padded_transcript` for the 3'UTR
/// side: a homopolymer that straddles the CDS↔3'UTR axis boundary so a
/// 5'-direction shuffle of `c.*N` would (without the axis clamp) shift
/// across into the CDS and re-classify as `c.<N>`.
///
/// Transcript layout:
///   tx 1-3   : 5'UTR `CCC`             (c.-3 .. c.-1)
///   tx 4-12  : CDS   `ATGCGCAAA`       (c.1 .. c.9; ends with AAA at c.7..c.9)
///   tx 13-20 : 3'UTR `AAAACCCC`        (c.*1 .. c.*8; starts with AAAA)
///
/// The combined A-tract runs tx 10-16 (c.7 .. c.*4) — 7 contiguous A's
/// straddling the CDS↔3'UTR boundary at tx 12/13 (c.9 / c.*1). Without
/// the clamp, a 5'-shuffle of `c.*1del` walks back through `AAAA` in the
/// 3'UTR and the `AAA` of the CDS, landing at `c.7del`. The clamp must
/// pin the 5'-direction shuffle at `c.*1` — the leftmost reachable
/// position on the 3'UTR axis.
fn provider_with_cds_utr3_homopolymer_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "CCCATGCGCAAAAAAACCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_UTR3.1".to_string(),
        Some("UTR3".to_string()),
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

fn normalize_utr3_with(direction: ShuffleDirection, cross: bool, input: &str) -> String {
    let mut config = NormalizeConfig::default().with_direction(direction);
    if cross {
        config = config.allow_crossing_boundaries();
    }
    let normalizer =
        Normalizer::with_config(provider_with_cds_utr3_homopolymer_transcript(), config);
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn five_prime_no_cross_does_not_shift_c_star_1_del_into_cds() {
    // 3'UTR mirror of `five_prime_no_cross_does_not_shift_c1_del_into_5utr`.
    // `c.*1del` under 5prime + cross_boundaries=false must stay at c.*1del
    // — the 5'-direction shuffle clamps at the CDS↔3'UTR boundary even
    // though the A-tract continues 3 bases into the CDS.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::FivePrime, false, "NM_UTR3.1:c.*1del"),
        "NM_UTR3.1:c.*1del",
    );
}

#[test]
fn five_prime_cross_does_not_shift_c_star_1_del_into_cds() {
    // Companion of `five_prime_cross_does_not_shift_c1_del_into_5utr`.
    // `cross_boundaries=true` authorizes crossing exon-intron junctions,
    // NOT the CDS↔3'UTR sub-axis boundary. The 3'UTR axis clamp fires
    // in both cross modes.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::FivePrime, true, "NM_UTR3.1:c.*1del"),
        "NM_UTR3.1:c.*1del",
    );
}

#[test]
fn five_prime_no_cross_does_not_shift_c_star_1_dup_into_cds() {
    // 3'UTR dup mirror of `five_prime_no_cross_does_not_shift_c1_dup_into_5utr`.
    // `c.*1dup` under 5prime + cross=false must stay at c.*1dup.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::FivePrime, false, "NM_UTR3.1:c.*1dup"),
        "NM_UTR3.1:c.*1dup",
    );
}

#[test]
fn five_prime_cross_does_not_shift_c_star_1_dup_into_cds() {
    // 3'UTR dup × cross=true cell.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::FivePrime, true, "NM_UTR3.1:c.*1dup"),
        "NM_UTR3.1:c.*1dup",
    );
}

#[test]
fn five_prime_inside_utr3_walks_to_c_star_1() {
    // `c.*4del` inside the 3'UTR's AAAA tract. 5'-direction shuffle
    // should walk through the 3'UTR A's and stop at `c.*1` — within
    // the 3'UTR axis, not crossing into CDS. Sanity check that the
    // clamp doesn't over-restrict (i.e. doesn't pin the variant at
    // its input position when the shuffle is entirely within-axis).
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::FivePrime, false, "NM_UTR3.1:c.*4del"),
        "NM_UTR3.1:c.*1del",
    );
}

#[test]
fn axis_clamp_emits_warning_when_constraining_5prime_c_star_1_shuffle() {
    // 3'UTR mirror of `axis_clamp_emits_warning_when_constraining_5prime_shuffle`.
    // Without the axis clamp `c.*1del` under 5prime would shift back
    // into the CDS A-tract (per fixture: CDS ends with `AAA`, 3'UTR
    // starts with `AAAA`). The clamp holds the result at c.*1del; the
    // clamp firing must emit `AxisClampApplied` so callers can flag.
    let normalizer = Normalizer::with_config(
        provider_with_cds_utr3_homopolymer_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_UTR3.1:c.*1del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_UTR3.1:c.*1del");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "AXIS_CLAMP_APPLIED"),
        "expected AXIS_CLAMP_APPLIED warning on 3'UTR 5'-shuffle that the \
         clamp held; got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
}

#[test]
fn axis_clamp_emits_warning_when_constraining_5prime_c_star_1_dup_shuffle() {
    // Dup mirror of `axis_clamp_emits_warning_when_constraining_5prime_c_star_1_shuffle`.
    // Same scenario, but with the dup edit kind to close the
    // (del × dup) leg of the warning matrix.
    let normalizer = Normalizer::with_config(
        provider_with_cds_utr3_homopolymer_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_UTR3.1:c.*1dup").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_UTR3.1:c.*1dup");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "AXIS_CLAMP_APPLIED"),
        "expected AXIS_CLAMP_APPLIED warning on 3'UTR 5'-shuffle of c.*1dup; \
         got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
}

// --- 3'-direction CDS-end clamp: c.<cds_end> 3'-shuffle must not leak --
//     into the 3'UTR. Mirror of the 5'-direction c.1 clamp on the
//     CDS-start side. The fixture re-uses `provider_with_cds_utr3_homopolymer_transcript`
//     because the CDS-end A-tract abuts the 3'UTR A-tract.

#[test]
fn three_prime_no_cross_does_not_shift_c_cds_end_del_into_utr3() {
    // 3'-direction symmetry of `five_prime_no_cross_does_not_shift_c1_del_into_5utr`.
    // `c.9del` (last CDS base, an A) under 3prime + cross_boundaries=false
    // must stay at c.9del — the 3'-direction shuffle clamps at the
    // CDS↔3'UTR boundary even though the A-tract continues 4 bases into
    // the 3'UTR.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::ThreePrime, false, "NM_UTR3.1:c.9del"),
        "NM_UTR3.1:c.9del",
    );
}

#[test]
fn three_prime_cross_does_not_shift_c_cds_end_del_into_utr3() {
    // `cross_boundaries=true` authorizes crossing exon-intron junctions
    // but NOT the CDS↔3'UTR sub-axis. Clamp fires in both cross modes.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::ThreePrime, true, "NM_UTR3.1:c.9del"),
        "NM_UTR3.1:c.9del",
    );
}

#[test]
fn three_prime_no_cross_does_not_shift_c_cds_end_dup_into_utr3() {
    // Dup symmetry of `three_prime_no_cross_does_not_shift_c_cds_end_del_into_utr3`.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::ThreePrime, false, "NM_UTR3.1:c.9dup"),
        "NM_UTR3.1:c.9dup",
    );
}

#[test]
fn three_prime_cross_does_not_shift_c_cds_end_dup_into_utr3() {
    // Closes the (del × dup) × (cross={false,true}) matrix on the
    // 3'-direction CDS-end side.
    assert_eq!(
        normalize_utr3_with(ShuffleDirection::ThreePrime, true, "NM_UTR3.1:c.9dup"),
        "NM_UTR3.1:c.9dup",
    );
}

#[test]
fn axis_clamp_emits_warning_when_constraining_3prime_c_cds_end_shuffle() {
    // 3'-direction mirror of the c.*1 5'-direction warning test. The
    // clamp prevents `c.9del` from shifting into the 3'UTR's AAAA
    // tract; that constraint must surface as `AXIS_CLAMP_APPLIED`.
    let normalizer = Normalizer::with_config(
        provider_with_cds_utr3_homopolymer_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs("NM_UTR3.1:c.9del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_UTR3.1:c.9del");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "AXIS_CLAMP_APPLIED"),
        "expected AXIS_CLAMP_APPLIED warning on CDS-end 3'-shuffle that the \
         clamp held; got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
}

/// Minus-strand 5'UTR-side fixture: identical tx-frame layout to
/// `provider_with_utr_padded_transcript` but on `Strand::Minus`.
///
/// The CDS↔UTR axis clamp operates on tx-frame positions
/// (`src/normalize/boundary.rs:get_cds_boundaries_with_axis_info`), so
/// the strand flag is by construction a no-op for the clamp itself.
/// These minus-strand tests therefore pin the **no-branch-on-strand**
/// contract — they're not exercising reverse-complement coordinate
/// arithmetic, just asserting that flipping `Strand::Plus` →
/// `Strand::Minus` does not introduce a stray strand-conditional code
/// path into the clamp.
fn provider_with_minus_strand_utr_padded_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "AAAATGAAATAGCCCCCCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_MTEST.1".to_string(),
        Some("MTEST".to_string()),
        Strand::Minus,
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

fn normalize_minus_with(direction: ShuffleDirection, cross: bool, input: &str) -> String {
    let mut config = NormalizeConfig::default().with_direction(direction);
    if cross {
        config = config.allow_crossing_boundaries();
    }
    let normalizer =
        Normalizer::with_config(provider_with_minus_strand_utr_padded_transcript(), config);
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn minus_strand_five_prime_no_cross_does_not_shift_c1_del_into_5utr() {
    // Minus-strand mirror of `five_prime_no_cross_does_not_shift_c1_del_into_5utr`.
    // The clamp must fire identically regardless of strand because it
    // operates on tx-frame positions.
    assert_eq!(
        normalize_minus_with(ShuffleDirection::FivePrime, false, "NM_MTEST.1:c.1del"),
        "NM_MTEST.1:c.1del",
    );
}

#[test]
fn minus_strand_five_prime_cross_does_not_shift_c1_del_into_5utr() {
    // Minus-strand × cross=true. Same strand-invariance contract.
    assert_eq!(
        normalize_minus_with(ShuffleDirection::FivePrime, true, "NM_MTEST.1:c.1del"),
        "NM_MTEST.1:c.1del",
    );
}

#[test]
fn minus_strand_five_prime_no_cross_does_not_shift_c1_dup_into_5utr() {
    // Minus-strand dup mirror of `five_prime_no_cross_does_not_shift_c1_dup_into_5utr`.
    assert_eq!(
        normalize_minus_with(ShuffleDirection::FivePrime, false, "NM_MTEST.1:c.1dup"),
        "NM_MTEST.1:c.1dup",
    );
}

#[test]
fn minus_strand_five_prime_cross_does_not_shift_c1_dup_into_5utr() {
    // Closes the minus-strand (del × dup) × (cross={false,true})
    // matrix on the 5'UTR side.
    assert_eq!(
        normalize_minus_with(ShuffleDirection::FivePrime, true, "NM_MTEST.1:c.1dup"),
        "NM_MTEST.1:c.1dup",
    );
}

#[test]
fn minus_strand_axis_clamp_emits_warning_when_constraining_5prime_shuffle() {
    // Minus-strand mirror of `axis_clamp_emits_warning_when_constraining_5prime_shuffle`.
    // The AXIS_CLAMP_APPLIED warning must fire on the minus-strand
    // fixture too — pinning that the warning-emission path is also
    // strand-invariant.
    let normalizer = Normalizer::with_config(
        provider_with_minus_strand_utr_padded_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_MTEST.1:c.1del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_MTEST.1:c.1del");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "AXIS_CLAMP_APPLIED"),
        "expected AXIS_CLAMP_APPLIED warning on minus-strand 5'-shuffle that \
         the clamp held; got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
}

#[test]
fn minus_strand_three_prime_shift_inside_cds_still_shuffles() {
    // Negative control: ensure the minus-strand path didn't accidentally
    // disable in-axis shuffling. `c.4del` (start of the CDS-internal
    // AAA tract) under 3'-direction must reach c.6del — the rightmost
    // CDS-internal A. Mirror of `three_prime_shift_inside_cds_still_shuffles`.
    assert_eq!(
        normalize_minus_with(ShuffleDirection::ThreePrime, false, "NM_MTEST.1:c.4del"),
        "NM_MTEST.1:c.6del",
    );
}

/// Minus-strand 3'UTR-side fixture: identical tx-frame layout to
/// `provider_with_cds_utr3_homopolymer_transcript` but on
/// `Strand::Minus`. Closes the minus-strand cell of the c.*1 / 3'UTR
/// matrix.
fn provider_with_minus_strand_cds_utr3_homopolymer_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "CCCATGCGCAAAAAAACCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_MUTR3.1".to_string(),
        Some("MUTR3".to_string()),
        Strand::Minus,
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

fn normalize_minus_utr3_with(direction: ShuffleDirection, cross: bool, input: &str) -> String {
    let mut config = NormalizeConfig::default().with_direction(direction);
    if cross {
        config = config.allow_crossing_boundaries();
    }
    let normalizer = Normalizer::with_config(
        provider_with_minus_strand_cds_utr3_homopolymer_transcript(),
        config,
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn minus_strand_five_prime_no_cross_does_not_shift_c_star_1_del_into_cds() {
    // Minus-strand mirror of `five_prime_no_cross_does_not_shift_c_star_1_del_into_cds`.
    assert_eq!(
        normalize_minus_utr3_with(ShuffleDirection::FivePrime, false, "NM_MUTR3.1:c.*1del"),
        "NM_MUTR3.1:c.*1del",
    );
}

#[test]
fn minus_strand_five_prime_cross_does_not_shift_c_star_1_del_into_cds() {
    // Closes minus-strand × 3'UTR × del × cross=true.
    assert_eq!(
        normalize_minus_utr3_with(ShuffleDirection::FivePrime, true, "NM_MUTR3.1:c.*1del"),
        "NM_MUTR3.1:c.*1del",
    );
}

#[test]
fn minus_strand_axis_clamp_emits_warning_when_constraining_5prime_c_star_1_shuffle() {
    // Minus-strand mirror of
    // `axis_clamp_emits_warning_when_constraining_5prime_c_star_1_shuffle`.
    let normalizer = Normalizer::with_config(
        provider_with_minus_strand_cds_utr3_homopolymer_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_MUTR3.1:c.*1del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("normalize");
    assert_eq!(format!("{}", result.result), "NM_MUTR3.1:c.*1del");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "AXIS_CLAMP_APPLIED"),
        "expected AXIS_CLAMP_APPLIED warning on minus-strand 3'UTR \
         5'-shuffle that the clamp held; got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
}
