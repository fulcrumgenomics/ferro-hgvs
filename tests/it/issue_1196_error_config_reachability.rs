//! Issue #1196 — error types that never reach a configuration consult site.
//!
//! Six `ErrorType` variants were never passed to `ErrorConfig::action_for` /
//! `should_reject` / `should_correct` / `should_warn`, and never attached to a
//! `DetectedCorrection` (which the preprocessor resolves through `action_for`).
//! For those codes `--error-mode` and the `--ignore` / `--reject` overrides
//! were accepted and then silently did nothing.
//!
//! The structural half of the fix lives in `error_code_audit.rs`, which now
//! scans for *consult* sites rather than for a bare mention of the variant
//! name. This file covers the two behavioral halves:
//!
//!   1. **W3022 `InitiatorMetCanonicalization` is now gated.** The advisory was
//!      pushed unconditionally, so silent mode and `--ignore W3022` could not
//!      suppress it. Every test below asserts its premise first, so none can
//!      decay into a tautology.
//!   2. **Genuinely non-configurable codes say so.** W3017/W3018/W3019 are hard
//!      parse rejections raised from the grammar, which runs below the
//!      error-configuration layer and has no correction to select; W2004/W3007/
//!      W3009 are registered but never emitted. All six now declare
//!      `consults_error_config() == false`, and the CLI reports an override
//!      naming one instead of accepting it and ignoring it.

use std::process::Command;

use ferro_hgvs::error_handling::{
    ErrorConfig, ErrorMode, ErrorOverride, ErrorType, ResolvedAction,
};
use ferro_hgvs::normalize::{NormalizationWarning, NormalizeConfig};
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// The `MKVAAALELE` fixture from #92's tests: inserting `Met` between p.1 and
/// p.2 canonicalizes to `p.Met1dup`, whose interval covers the initiator
/// methionine, so W3022 fires.
fn protein_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein("NP_TESTPROT.1", "MKVAAALELE");
    provider
}

/// Normalize `input` under `error_config` and report whether the W3022
/// advisory was surfaced.
fn met1_advisory_emitted(input: &str, error_config: ErrorConfig) -> bool {
    let config = NormalizeConfig::default().with_error_config(error_config);
    let normalizer = Normalizer::with_config(protein_provider(), config);
    let variant = parse_hgvs(input).expect("fixture input must parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("W3022 is advisory — normalization must succeed in every mode");
    result
        .warnings
        .iter()
        .any(|w| matches!(w, NormalizationWarning::InitiatorMetCanonicalization { .. }))
}

const MET1_INS: &str = "NP_TESTPROT.1:p.Met1_Lys2insMet";

// ---------------------------------------------------------------------------
// W3022 — the advisory now honors the error configuration
// ---------------------------------------------------------------------------

/// Strict and lenient both surface the advisory. Strict must not be *quieter*
/// than lenient, which is what a naive `should_warn()`-only gate would produce
/// (`ResolvedAction::should_warn` is true only for `WarnCorrect`, i.e. lenient).
#[test]
fn met1_advisory_surfaced_in_strict_and_lenient() {
    assert!(
        met1_advisory_emitted(MET1_INS, ErrorConfig::strict()),
        "strict mode must surface W3022",
    );
    assert!(
        met1_advisory_emitted(MET1_INS, ErrorConfig::lenient()),
        "lenient mode must surface W3022",
    );
}

/// The regression: silent mode must suppress the advisory. Before #1196 the
/// push was unconditional, so this returned `true` in every mode.
#[test]
fn met1_advisory_suppressed_in_silent_mode() {
    assert!(
        met1_advisory_emitted(MET1_INS, ErrorConfig::lenient()),
        "premise: the advisory fires for this input under lenient mode",
    );
    assert!(
        !met1_advisory_emitted(MET1_INS, ErrorConfig::silent()),
        "silent mode must suppress W3022",
    );
}

/// The `--ignore` direction: an `ErrorOverride::Accept` must suppress the
/// advisory even from a strict base, proving the override is read rather than
/// the mode being hardwired.
#[test]
fn met1_advisory_suppressed_by_an_ignore_override() {
    assert!(
        met1_advisory_emitted(MET1_INS, ErrorConfig::strict()),
        "premise: strict surfaces the advisory without an override",
    );
    let ignored = ErrorConfig::strict().with_override(
        ErrorType::InitiatorMetCanonicalization,
        ErrorOverride::Accept,
    );
    assert!(
        !met1_advisory_emitted(MET1_INS, ignored),
        "an Accept override (`--ignore W3022`) must suppress the advisory",
    );
}

/// The `--reject` direction must not be inert either.
///
/// An explicit `Reject` override — and *only* an explicit one — promotes the
/// advisory to a hard error. Without this, `--reject W3022` would resolve to the
/// same `ResolvedAction::Reject` as strict mode, be treated as "surface it", and
/// silently do nothing: exactly the defect #1196 exists to remove, reintroduced
/// for one code.
#[test]
fn met1_advisory_rejects_under_an_explicit_reject_override() {
    let config =
        NormalizeConfig::default().with_error_config(ErrorConfig::lenient().with_override(
            ErrorType::InitiatorMetCanonicalization,
            ErrorOverride::Reject,
        ));
    assert!(
        config.should_reject_initiator_met_canonicalization(),
        "an explicit Reject override must be visible to the normalizer",
    );
    assert!(
        config.should_warn_initiator_met_canonicalization(),
        "premise: the advisory is still surfaced, because the rejection ladder promotes it \
         from the warning list",
    );

    let normalizer = Normalizer::with_config(protein_provider(), config);
    let variant = parse_hgvs(MET1_INS).expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("`--reject W3022` must reject");
    let rendered = format!("{err}");
    assert!(
        rendered.contains("W3022"),
        "the rejection must name the code it was asked to reject; got {rendered:?}",
    );
}

/// The base mode alone must NEVER promote the advisory — the distinction between
/// "strict" and "explicitly rejected" is what keeps strict mode from refusing
/// ferro's own canonical output.
#[test]
fn met1_advisory_is_not_promoted_by_strict_mode_alone() {
    let strict = NormalizeConfig::default().with_error_config(ErrorConfig::strict());
    assert_eq!(
        strict.initiator_met_canonicalization_action(),
        ResolvedAction::Reject,
        "premise: strict resolves to Reject, so the gate cannot key on the resolved action \
         alone",
    );
    assert!(
        !strict.should_reject_initiator_met_canonicalization(),
        "strict mode must not promote an advisory about ferro's own output",
    );
    assert!(
        strict.should_warn_initiator_met_canonicalization(),
        "strict must instead surface it, so it is never quieter than lenient",
    );
}

/// W3022 must never become a rejection, in any mode.
///
/// The warning fires whenever the *final* edit is a Met1-covering duplication —
/// including when the input already was one and nothing was rewritten. It is an
/// advisory about ferro's own canonical output, so promoting it to an error (as
/// the registry's former `warn_accept` mode table implied) would make strict
/// mode refuse the very string strict mode emits, and would break
/// normalization idempotency. This pins the decision.
#[test]
fn met1_advisory_never_rejects_ferros_own_canonical_output() {
    let canonical = "NP_TESTPROT.1:p.Met1dup";
    for (label, error_config) in [
        ("strict", ErrorConfig::strict()),
        ("lenient", ErrorConfig::lenient()),
        ("silent", ErrorConfig::silent()),
    ] {
        let config = NormalizeConfig::default().with_error_config(error_config);
        let normalizer = Normalizer::with_config(protein_provider(), config);
        let variant = parse_hgvs(canonical).expect("parse p.Met1dup");
        let normalized = normalizer
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("{label} mode must accept ferro's own output: {e}"));
        assert_eq!(
            format!("{normalized}"),
            canonical,
            "{label} mode must round-trip the canonical form unchanged",
        );
    }
}

// ---------------------------------------------------------------------------
// The declaration, and the CLI diagnostic it drives
// ---------------------------------------------------------------------------

/// The codes whose handling never consults the error configuration.
///
/// Three disjoint reasons, pinned together because the user-visible consequence
/// is the same: an override naming any of them cannot change anything. The
/// count is deliberately not stated — it has already gone stale once, and a
/// number in a doc comment guards nothing.
#[test]
fn non_configurable_codes_are_declared_inert() {
    // Registered but never emitted — retained as stable Python discriminants
    // (#1114).
    for error_type in [
        ErrorType::InvalidUnicodeCharacter,
        ErrorType::DeprecatedStopCodonStar,
        ErrorType::DeprecatedFrameshiftStar,
    ] {
        assert!(
            !error_type.consults_error_config(),
            "{error_type:?} is never emitted, so no mode can act on it",
        );
    }

    // Hard parse rejections raised below the error-configuration layer, with no
    // auto-correction to select.
    for error_type in [
        ErrorType::AlleleFractionAnnotation,
        ErrorType::ClinVarProseMultiAllelic,
        ErrorType::NonSpecMosaicForm,
    ] {
        assert!(
            !error_type.consults_error_config(),
            "{error_type:?} rejects in every mode",
        );
        assert!(
            !error_type.is_correctable(),
            "premise: {error_type:?} has no auto-correction, which is why no mode \
             can soften it",
        );
    }

    // Emitted unconditionally by the normalizer, with no consult site: the
    // warning reports provenance about ferro's own output rather than a defect
    // in the input, so there is nothing to reject and nothing to correct
    // (#2092). Its registry entry carries no `ModeBehavior` for the same
    // reason, which `registry::tests::test_warning_codes_have_mode_behavior`
    // keys on — so the two statements cannot drift apart silently.
    assert!(
        !ErrorType::MembersCoalescedFromReportedForm.consults_error_config(),
        "W5005 is pushed at the normalization exit in every mode, so no override reaches it",
    );
    assert!(
        !ErrorType::MembersCoalescedFromReportedForm.is_correctable(),
        "premise: the coalesced form IS the correct description, so there is nothing to correct",
    );
}

/// A representative configurable code must still declare itself reachable, so
/// the test above is not passing merely because the method returns `false`
/// everywhere.
#[test]
fn configurable_codes_are_declared_reachable() {
    for error_type in [
        ErrorType::RefSeqMismatch,
        ErrorType::DeprecatedStopCodonX,
        ErrorType::InitiatorMetCanonicalization,
    ] {
        assert!(
            error_type.consults_error_config(),
            "{error_type:?} is consulted at runtime and must declare it",
        );
    }
}

/// `ferro parse --ignore W3017` used to be accepted in silence: the code
/// resolved, the override was stored, and nothing ever read it. It must now say
/// so — while still rejecting the input, because the override genuinely cannot
/// change the outcome.
#[test]
fn cli_reports_an_inert_ignore_override() {
    let output = ferro()
        .args(["parse", "NC_012920.1:m.1234A>G(80%)", "--ignore", "W3017"])
        .output()
        .expect("run ferro parse");
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        stderr.contains("W3017"),
        "the inert-override diagnostic must name the code; stderr: {stderr}",
    );
    assert!(
        stderr.contains("has no effect"),
        "the diagnostic must say the override does nothing; stderr: {stderr}",
    );
    assert!(
        stderr.contains("Allele-fraction"),
        "the input must still be rejected by the parser; stderr: {stderr}",
    );
}

/// The companion premise: a *configurable* code must NOT draw the inert
/// diagnostic, so the warning cannot degrade into unconditional noise on every
/// `--ignore`.
#[test]
fn cli_does_not_report_a_working_override_as_inert() {
    let output = ferro()
        .args([
            "parse",
            "NM_000088.3:c.100A>G",
            "--error-mode",
            "strict",
            "--ignore",
            "W2001",
        ])
        .output()
        .expect("run ferro parse");
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        !stderr.contains("has no effect"),
        "W2001 (WrongDashCharacter) is consulted by the preprocessor and must not be \
         reported as inert; stderr: {stderr}",
    );
}

/// Every `ErrorMode` resolves an inert code the same way, which is the property
/// that makes the CLI diagnostic true rather than merely cautious.
#[test]
fn inert_codes_reject_under_every_mode() {
    let input = "NC_012920.1:m.1234A>G(80%)";
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let config = ErrorConfig::new(mode)
            .with_override(ErrorType::AlleleFractionAnnotation, ErrorOverride::Accept);
        let result = ferro_hgvs::hgvs::parser::parse_hgvs_with_config(input, config);
        assert!(
            result.is_err(),
            "{mode:?} + an Accept override must still reject {input}: the code is \
             declared inert precisely because nothing reads that override",
        );
    }
}
