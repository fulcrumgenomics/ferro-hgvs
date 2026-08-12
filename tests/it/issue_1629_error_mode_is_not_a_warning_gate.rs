//! #1629 — `ErrorMode` carries no mode-level warning gate, and must not acquire one.
//!
//! `ErrorMode::emits_warnings()` was a `matches!(self, Lenient)` that shipped in
//! the initial commit and never had a call site anywhere in the crate. `git log
//! -S emits_warnings -- src/` returns exactly one commit — that initial one — so
//! this was never a caller that got removed; it was a contract stated and never
//! enforced. Its sibling `allows_correction()` was the same shape and equally
//! unreferenced. Both are gone.
//!
//! # The tempting repair, and why it is wrong
//!
//! #1629 offered two shapes: delete it, or "give `emits_warnings()` a call site
//! — one gate at the point warnings are collected, so the mode decides and no
//! code can opt out by omission". **The second is not available**, and the tests
//! below are what say so rather than a paragraph asserting it.
//!
//! Warning emission is resolved *per code*, by
//! `ErrorConfig::action_for(code) = overrides[code].resolve(mode)`, and read
//! through `ResolvedAction::should_warn`. The mode is only the base of that
//! resolution. `--ignore <code>` and `--reject <code>` (#1196) exist precisely
//! to move one code off its mode's default, in **both** directions:
//!
//! - a `WarnCorrect` override under `Silent` warns — the operator asked for it;
//! - an `Accept`/`SilentCorrect` override under `Lenient` does not.
//!
//! A gate keyed on the mode would overrule the first and be redundant with the
//! second, i.e. it would delete the knob it was placed next to. That is a worse
//! defect than the dead predicate.
//!
//! Ferro also has a code whose emission deliberately does **not** follow
//! `Lenient`: `InitiatorMetCanonicalization` (W3022) is an advisory about
//! ferro's *own* canonical output, so promoting it to an error in strict mode
//! would make strict refuse the string strict emits. Its resolution is therefore
//! `should_warn() || should_reject()`, which is true under `Strict` — where a
//! `Strict.emits_warnings() == false` gate would have suppressed it. The
//! deleted predicate was not merely unused; on this code it was wrong.
//!
//! # What survives
//!
//! `ErrorMode::is_strict()` does, because it is a different kind of question: it
//! reports the mode, and its one caller
//! (`NormalizeConfig::should_reject_reduced_capability`) is documented as
//! deliberately *not* on the per-code errors axis — `ReducedCapabilityNoGenome`
//! is an environmental limitation with no `ErrorType` for an override to name.
//! It is pinned below against overrides for exactly that reason.
//!
//! # Not covered here, and it is an operator decision
//!
//! Two per-code predicates are themselves unwired —
//! `NormalizeConfig::should_warn_ref_mismatch` and
//! `should_warn_variant_exceeds_reference` have no call site outside their own
//! unit tests, and `src/normalize/mod.rs` pushes `RefSeqMismatch`
//! unconditionally. Wiring them would make `--error-mode silent` suppress
//! `W5001` (`RefSeqMismatch`) and `W5003` (`VariantExceedsReference`) — the two
//! codes those predicates actually resolve; `W3001` is `MissingVersion`, an
//! unrelated preprocessor code. That is a behaviour change on a question
//! `issue_1181_cli_error_mode.rs::strict_is_distinct_from_lenient_and_silent_at_the_cli`
//! records as **open** (it pins `lenient_out == silent_out` today). Deciding it
//! is an adjudication, not a cleanup, so nothing here touches it.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode, ErrorOverride, ErrorType};
use ferro_hgvs::NormalizeConfig;

/// Every code that could plausibly be gated by a mode-level warning predicate,
/// so the properties below are not statements about one lucky variant.
const OVERRIDABLE_CODES: &[ErrorType] = &[
    ErrorType::WrongDashCharacter,
    ErrorType::LowercaseAminoAcid,
    ErrorType::MissingVersion,
    ErrorType::ExtraWhitespace,
    ErrorType::TrailingAnnotation,
];

/// **Question.** Does the mode alone decide whether a code warns?
///
/// **No, in both directions.** This is the property that makes a mode-level
/// `emits_warnings()` unwireable rather than merely unused: for every
/// overridable code, `Silent` can be made to warn and `Lenient` can be made not
/// to. A gate on the mode would have to ignore one of these to exist.
#[test]
fn an_override_moves_warning_emission_off_the_mode_in_both_directions() {
    for &code in OVERRIDABLE_CODES {
        let silent_default = ErrorConfig::silent();
        assert!(
            !silent_default.should_warn(code),
            "{code:?}: silent mode's default must not warn — the premise of the two flips below"
        );
        let silent_but_warning =
            ErrorConfig::silent().with_override(code, ErrorOverride::WarnCorrect);
        assert!(
            silent_but_warning.should_warn(code),
            "{code:?}: `--reject`/`--warn`-style override must reach a code under SILENT mode; a \
             mode-level warning gate would silence it and delete the knob"
        );

        let lenient_default = ErrorConfig::lenient();
        assert!(
            lenient_default.should_warn(code),
            "{code:?}: lenient mode's default must warn — the premise of the flip below"
        );
        let lenient_but_silent =
            ErrorConfig::lenient().with_override(code, ErrorOverride::SilentCorrect);
        assert!(
            !lenient_but_silent.should_warn(code),
            "{code:?}: `--ignore <code>` must silence one code under LENIENT mode"
        );
    }
}

/// **Question.** Is `Lenient` the only mode that emits a warning, as the deleted
/// predicate said?
///
/// **No.** W3022 is emitted under `Strict`, deliberately and with a recorded
/// reason (`NormalizeConfig::should_warn_initiator_met_canonicalization`): strict
/// mode must not refuse ferro's own canonical output, so `Reject` from the base
/// mode maps to "surface it" rather than to "error". `Strict.emits_warnings()`
/// answered `false`, so the gate #1629 proposed would have suppressed a
/// diagnostic that exists to keep strict mode from being *quieter* than lenient.
#[test]
fn strict_mode_emits_the_initiator_met_advisory() {
    assert!(
        NormalizeConfig::strict().should_warn_initiator_met_canonicalization(),
        "W3022 must surface under strict mode; a `Strict => no warnings` gate would suppress it"
    );
    assert!(
        NormalizeConfig::lenient().should_warn_initiator_met_canonicalization(),
        "…and under lenient"
    );
    assert!(
        !NormalizeConfig::silent().should_warn_initiator_met_canonicalization(),
        "…and silent is where it is actually suppressed, per-code, by the resolution"
    );
}

/// **Question.** Does the one surviving mode-level predicate stay a statement
/// about the mode?
///
/// **Yes, and that is what makes it wireable where the other two were not.**
/// `should_reject_reduced_capability` asks a question with no `ErrorType`
/// attached — `ReducedCapabilityNoGenome` is an environmental limitation, not an
/// input defect (#1012 item 2) — so no override can or should move it. Pinned
/// against a config carrying overrides so that re-expressing it through
/// `action_for(..)` (the change this test exists to block) fails here.
#[test]
fn is_strict_is_wired_and_answers_the_mode_not_a_code() {
    let strict_with_overrides = OVERRIDABLE_CODES
        .iter()
        .fold(ErrorConfig::strict(), |config, &code| {
            config.with_override(code, ErrorOverride::Accept)
        });
    assert!(
        NormalizeConfig::default()
            .with_error_config(strict_with_overrides)
            .should_reject_reduced_capability(),
        "reduced capability is off the per-code axis: no `--ignore` may switch it off"
    );

    for mode in [ErrorMode::Lenient, ErrorMode::Silent] {
        let rejecting_everything = OVERRIDABLE_CODES
            .iter()
            .fold(ErrorConfig::new(mode), |config, &code| {
                config.with_override(code, ErrorOverride::Reject)
            });
        assert!(
            !NormalizeConfig::default()
                .with_error_config(rejecting_everything)
                .should_reject_reduced_capability(),
            "{mode}: nor may a `--reject` switch it on"
        );
        assert!(!mode.is_strict(), "{mode} is not strict");
    }
    assert!(ErrorMode::Strict.is_strict());
}
