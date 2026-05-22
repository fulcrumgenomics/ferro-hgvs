//! Issue #330 — structured info-code surface mirroring the W#### warning
//! surface.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/330>.
//!
//! Production users currently cannot distinguish "ferro silently shifted
//! your variant per the 3' rule" from "ferro left it alone" without
//! re-comparing strings. This pins the canonical contract:
//!
//! - `Normalizer::normalize_with_diagnostics` returns a result carrying both
//!   `warnings: Vec<NormalizationWarning>` and `infos: Vec<NormalizationInfo>`.
//! - Each [`NormalizationInfo`] exposes `code()` (a stable identifier such as
//!   `"SHUFFLE_APPLIED"`) and `message()` (a human-readable description).
//! - When the shuffle layer relocates a variant, `infos` contains a
//!   [`NormalizationInfo::ShuffleApplied`] variant carrying the direction,
//!   the original HGVS position text, and the normalized HGVS position
//!   text.
//! - When the variant is already canonical (or the edit kind is not
//!   shuffle-eligible), `infos` is empty.
//!
//! Tests use a `MockProvider`-backed `c.` variant on a homopolymer tract,
//! mirroring the example given in the issue body. Complementary tests pin
//! the no-info empty case, multi-position spans, the 5'-shift direction,
//! and that the non-`with_warnings` `normalize()` entry point
//! intentionally discards infos. The conservative no-op behaviour on
//! top-level kind changes (`HV::Allele` ↔ bare variant across
//! cis-collapse / `wrap_allele_if_split`) is documented on
//! `detect_shuffle_infos` itself rather than pinned here.

use ferro_hgvs::normalize::config::ShuffleDirection;
use ferro_hgvs::normalize::{NormalizationInfo, NormalizeConfig, Normalizer};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};

/// Build a single-exon plus-strand transcript whose CDS contains a run of
/// `A`s, so that `c.<pos>delA` inside the run is shuffle-eligible.
fn provider_with_polya() -> MockProvider {
    let mut provider = MockProvider::new();
    // Position:  1234567890123456789012345
    // Sequence:  ATGAAAAAAAAACGTCGTCGTCGTC
    //               ^^^^^^^^^^ poly-A tract (c.4..c.12)
    provider.add_transcript(Transcript::new(
        "NM_POLYA.1".to_string(),
        Some("POLYA".to_string()),
        Strand::Plus,
        Some("ATGAAAAAAAAACGTCGTCGTCGTC".to_string()),
        Some(1),
        Some(25),
        vec![Exon::new(1, 1, 25)],
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

/// `c.4del` (a single-base deletion inside the poly-A tract) must 3'-shift
/// to `c.12del` (the rightmost canonical position). The diagnostics surface
/// records the shift via a single [`NormalizationInfo::ShuffleApplied`]
/// with `direction = ThreePrime`.
#[test]
fn shuffle_emits_info_on_homopolymer_three_prime() {
    let normalizer = Normalizer::new(provider_with_polya());
    let variant = parse_hgvs("NM_POLYA.1:c.4del").expect("parse c.4del");

    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize c.4del");

    let display = format!("{}", result.result);
    assert!(
        display.ends_with(":c.12del"),
        "expected shuffle to land on c.12del; got {display:?}",
    );

    assert_eq!(
        result.infos.len(),
        1,
        "expected one info code for the shift; got {:?}",
        result.infos.iter().map(|i| i.code()).collect::<Vec<_>>(),
    );
    let info = &result.infos[0];
    assert_eq!(info.code(), "SHUFFLE_APPLIED");
    let NormalizationInfo::ShuffleApplied {
        direction,
        original_position,
        normalized_position,
        ..
    } = info
    else {
        panic!("expected ShuffleApplied, got {info:?}");
    };
    assert_eq!(*direction, ShuffleDirection::ThreePrime);
    assert_eq!(original_position, "4");
    assert_eq!(normalized_position, "12");
}

/// 5'-shift configuration: a deletion at the right edge of the poly-A
/// tract is shuffled leftward to its 5'-most position. The emitted info
/// must carry `direction = FivePrime` so callers can interpret the
/// signal correctly.
#[test]
fn shuffle_info_carries_direction_under_five_prime_config() {
    let config = NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime);
    let normalizer = Normalizer::with_config(provider_with_polya(), config);
    let variant = parse_hgvs("NM_POLYA.1:c.12del").expect("parse c.12del");

    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize c.12del under 5'-shift");

    let display = format!("{}", result.result);
    assert!(
        display.ends_with(":c.4del"),
        "expected 5'-shuffle to land on c.4del; got {display:?}",
    );

    assert_eq!(result.infos.len(), 1, "one info expected");
    let NormalizationInfo::ShuffleApplied {
        direction,
        original_position,
        normalized_position,
        ..
    } = &result.infos[0]
    else {
        panic!("expected ShuffleApplied, got {:?}", result.infos[0]);
    };
    assert_eq!(
        *direction,
        ShuffleDirection::FivePrime,
        "info must record the configured direction, not assume 3'",
    );
    assert_eq!(original_position, "12");
    assert_eq!(normalized_position, "4");
}

/// Multi-position spans (`c.4_5del` inside a poly-A run) shuffle as a
/// unit. The diagnostics surface records the span move, with
/// `original_position` and `normalized_position` carrying the multi-base
/// `start_end` form.
#[test]
fn shuffle_info_for_multi_position_span() {
    let normalizer = Normalizer::new(provider_with_polya());
    let variant = parse_hgvs("NM_POLYA.1:c.4_5del").expect("parse c.4_5del");

    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize c.4_5del");

    let display = format!("{}", result.result);
    assert!(
        display.ends_with(":c.11_12del"),
        "expected 3'-shuffle to land on c.11_12del; got {display:?}",
    );

    assert_eq!(result.infos.len(), 1);
    let NormalizationInfo::ShuffleApplied {
        original_position,
        normalized_position,
        ..
    } = &result.infos[0]
    else {
        panic!("expected ShuffleApplied, got {:?}", result.infos[0]);
    };
    assert_eq!(original_position, "4_5");
    assert_eq!(normalized_position, "11_12");
}

/// `c.3G>C` on a unique position cannot shuffle (substitutions don't
/// shift, and the position is unique anyway). The diagnostics surface
/// must report an empty `infos` list.
#[test]
fn no_info_when_variant_is_already_canonical() {
    let normalizer = Normalizer::new(provider_with_polya());
    let variant = parse_hgvs("NM_POLYA.1:c.3G>C").expect("parse c.3G>C");

    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize c.3G>C");

    assert!(
        result.infos.is_empty(),
        "expected empty infos for canonical input; got {:?}",
        result.infos.iter().map(|i| i.code()).collect::<Vec<_>>(),
    );
}

/// `Normalizer::normalize()` (the non-`with_warnings` entry point)
/// returns only the normalized variant; infos are intentionally
/// discarded. Callers that need the diagnostic surface must use
/// `normalize_with_diagnostics`. Pin this contract so a future refactor
/// doesn't accidentally surface infos through the bare `normalize()`
/// return type.
#[test]
fn normalize_without_warnings_discards_infos_by_design() {
    let normalizer = Normalizer::new(provider_with_polya());
    let variant = parse_hgvs("NM_POLYA.1:c.4del").expect("parse c.4del");

    // Bare normalize returns just HgvsVariant — no infos surface.
    let normalized = normalizer.normalize(&variant).expect("normalize c.4del");
    // Sanity: the shuffle still happened.
    let display = format!("{}", normalized);
    assert!(
        display.ends_with(":c.12del"),
        "shuffle should still apply via bare normalize; got {display:?}",
    );

    // And the with_warnings entry point is the only path that surfaces
    // the info.
    let detailed = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize_with_diagnostics c.4del");
    assert_eq!(detailed.infos.len(), 1);
}

/// `NormalizationInfo::code()` returns a stable identifier suitable for
/// programmatic dispatch. Pin the spelling so downstream callers can match
/// on it without depending on `Debug` repr.
#[test]
fn info_code_is_stable_identifier() {
    let info = NormalizationInfo::ShuffleApplied {
        accession: "NM_POLYA.1".to_string(),
        direction: ShuffleDirection::ThreePrime,
        original_position: "4".to_string(),
        normalized_position: "12".to_string(),
    };
    assert_eq!(info.code(), "SHUFFLE_APPLIED");
    // `message()` synthesizes from structural fields (#397 item 3
    // dropped the per-variant `message: String` field).
    let msg = info.message();
    assert!(
        msg.contains("4") && msg.contains("12") && msg.contains("NM_POLYA.1"),
        "synthesized message should mention positions and accession; got {msg:?}",
    );
}
