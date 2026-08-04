//! #1382 — both public normalization exits must apply the same strict-mode
//! rejection ladder.
//!
//! `Normalizer` has exactly two public normalizing methods, and they routed
//! differently:
//!
//! | method | route | strict ladder |
//! | --- | --- | --- |
//! | `normalize()` | `normalize_core_checked` | yes |
//! | `normalize_with_diagnostics()` | `normalize_core_canonical` | **no** |
//!
//! `normalize_core_checked` is where the `should_reject_*` ladder lives —
//! `RefSeqMismatch`, W5002/W5003/W5004/W4004/W4005/W4006, `EINTRONIC`, and the
//! `ReducedCapabilityNoGenome` promotion. `normalize_core_canonical` is just
//! `normalize_core` plus `canonicalize_from_sequence`, so the diagnostics exit
//! returned a result no rejection had ever seen. A strict-configured
//! `Normalizer` therefore rejected a variant through one entry point and
//! accepted it through the other.
//!
//! It stayed hidden because in the default lenient config every
//! `should_reject_*` is false and the ladder is a no-op — the two exits agree
//! everywhere except under a strict config, and almost nothing normalizes
//! strictly through the diagnostics exit.
//!
//! ## How it surfaced
//!
//! `issue_336_position_past_end::strict_mode_n_in_bounds_does_not_emit_warning`
//! was passing a variant whose stated reference base is wrong — fixture
//! `ACGTACGTACGTACGTACGT`, where `n.10` is `C`, against an authored
//! `NR_TEST.1:n.10G>C` — under `NormalizeConfig::strict()`. The `RefSeqMismatch`
//! was never promoted on that path. #1366 corrected that fixture (the test is
//! about the position, not the base) and left the asymmetry for this issue.
//!
//! ## What these tests pin
//!
//! Agreement between the two exits, asserted as agreement rather than as two
//! independent expectations: each case runs the *same* variant through both and
//! requires the verdicts to match. A test that asserted "both reject" would go
//! green if a future change made `normalize()` stop rejecting, which is the
//! opposite of what this issue is about.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

/// A 20-base non-coding transcript, `ACGTACGT…`, so 1-based position `p` holds
/// `"ACGT"[(p - 1) % 4]`. Position 10 is `C` and position 20 is `T`.
///
/// Non-coding (no `cds_start`/`cds_end`) so the `n.` axis is the transcript's
/// own and no CDS arithmetic enters the picture.
fn noncoding_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "ACGTACGTACGTACGTACGT".to_string();
    let len = sequence.len() as u64;
    provider.add_transcript(Transcript::new(
        "NR_TEST.1".to_string(),
        Some("NCRNA".to_string()),
        Strand::Plus,
        sequence,
        None,
        None,
        vec![Exon::new(1, 1, len)],
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

/// Run one descriptor through both public exits under `config` and return
/// whether each rejected it.
fn verdicts(descriptor: &str, config: NormalizeConfig) -> (bool, bool) {
    let variant = parse_hgvs(descriptor).expect("fixture must parse");
    let normalizer = Normalizer::with_config(noncoding_provider(), config);
    (
        normalizer.normalize(&variant).is_err(),
        normalizer.normalize_with_diagnostics(&variant).is_err(),
    )
}

/// The reproducer. A wrong stated reference base under a strict config is a
/// `RefSeqMismatch`, which `normalize()` promotes to a hard error — so the
/// diagnostics exit must promote it too.
#[test]
fn a_strict_config_rejects_a_reference_mismatch_through_both_exits() {
    let (plain, diagnostics) = verdicts("NR_TEST.1:n.10G>C", NormalizeConfig::strict());
    assert!(
        plain,
        "`normalize()` must reject a stated `G` where the transcript has `C` under a \
         strict config — if this is false the premise of #1382 has changed and the \
         assertion below is measuring nothing"
    );
    assert_eq!(
        plain, diagnostics,
        "the two public exits disagree: `normalize()` rejected `n.10G>C` and \
         `normalize_with_diagnostics()` did not. Both must apply the same \
         strict-mode rejection ladder (#1382)"
    );
}

/// The same input past the transcript's end, so a second rung of the ladder is
/// covered rather than only `RefSeqMismatch`. `n.21` does not exist on a
/// 20-base transcript (W4004 `PositionPastEnd`).
#[test]
fn a_strict_config_rejects_a_past_end_position_through_both_exits() {
    let (plain, diagnostics) = verdicts("NR_TEST.1:n.21G>C", NormalizeConfig::strict());
    assert!(
        plain,
        "`normalize()` must reject `n.21` on a 20-base transcript under a strict config"
    );
    assert_eq!(
        plain, diagnostics,
        "the two public exits disagree on a past-end position: `normalize()` \
         rejected `n.21G>C` and `normalize_with_diagnostics()` did not (#1382)"
    );
}

/// Lenient is the default and must stay permissive on **both** exits. This is
/// the control: the fix routes the ladder through the diagnostics exit, and if
/// it accidentally rejected in lenient mode it would break every caller.
#[test]
fn the_default_lenient_config_still_accepts_on_both_exits() {
    for descriptor in ["NR_TEST.1:n.10G>C", "NR_TEST.1:n.10C>G"] {
        let (plain, diagnostics) = verdicts(descriptor, NormalizeConfig::default());
        assert!(
            !plain && !diagnostics,
            "`{descriptor}` must be accepted by both exits in the default lenient \
             config; got normalize()={plain:?} normalize_with_diagnostics()={diagnostics:?}"
        );
    }
}

/// A correct description must be accepted under a strict config through both
/// exits, so the fix is not simply making the diagnostics exit reject more.
#[test]
fn a_strict_config_accepts_a_correct_description_through_both_exits() {
    let (plain, diagnostics) = verdicts("NR_TEST.1:n.10C>G", NormalizeConfig::strict());
    assert!(
        !plain && !diagnostics,
        "`n.10C>G` states the real base at position 10 and is in bounds, so a strict \
         config must accept it on both exits; got normalize()={plain:?} \
         normalize_with_diagnostics()={diagnostics:?}"
    );
}

/// The warnings the diagnostics exit exists to surface must survive the change.
///
/// Routing through the ladder must not swallow the advisory warnings a lenient
/// caller reads — that is the whole point of this entry point, and a fix that
/// returned the variant with an empty warning list would pass every assertion
/// above while destroying the method's purpose.
#[test]
fn the_diagnostics_exit_still_reports_its_warnings() {
    let variant = parse_hgvs("NR_TEST.1:n.10G>C").expect("fixture must parse");
    let normalizer = Normalizer::with_config(
        noncoding_provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Lenient),
    );
    let diagnosed = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient must not reject");
    // The specific warning, not merely "some warning": this exit now runs the
    // whole ladder, which can surface unrelated warnings, and any one of them
    // would satisfy a non-empty check while the reference-mismatch report this
    // test exists for had quietly stopped being emitted.
    assert!(
        diagnosed
            .warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. })),
        "a stated `G` where the transcript has `C` must surface a RefSeqMismatch on \
         the diagnostics exit; got {:?}",
        diagnosed.warnings
    );
}
