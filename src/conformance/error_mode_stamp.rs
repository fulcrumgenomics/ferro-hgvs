//! The error mode an artifact was generated under, stamped into the artifact (#1629).
//!
//! # The defect
//!
//! Every census and corpus in this repository is a measurement taken through a
//! normalizer, and a normalizer carries an [`ErrorConfig`]. The mode that config
//! names is not a detail of the run — it changes the numbers. Measured on the
//! spec conformance corpus (`tests/it/corpus_prohibited_inputs.rs`), which
//! re-measures three of the four refusal counters in strict mode: two of them
//! move — `conflicts_accepted`'s 72 are all refused and 16 of
//! `prohibited_conditional_accepted`'s 40 are — while `prohibited_absolute_accepted`'s
//! 32 stay accepted. One counter that moves is already enough to make a census
//! compared across modes a category error, and the decided ruling
//! `bare-transcript-intronic-position` has to end with the sentence "**so that
//! counter is a lenient-mode figure**" because nothing in the artifact says so.
//!
//! The failure is quiet in the way that costs the most. `NormalizeConfig::default()`
//! substitutes `ErrorConfig::lenient()`, while `ErrorMode`'s own `#[default]` is
//! `Strict` and `ErrorConfig::default()` is `strict()` — so a generator written
//! as `Normalizer::new(provider)` measures in **lenient** mode while reading as
//! though it took the default, which everywhere else in the crate means strict.
//! `NormalizeConfig` cannot help: its `error_config` field is `#[serde(skip)]`,
//! so serializing the whole config still emits no mode.
//!
//! Two censuses taken under different modes are then indistinguishable
//! artifacts, and comparing them is a category error nobody can detect from the
//! files. That is the same shape as the corpus-zero hazard
//! [`completeness`](super::completeness) exists for: a measurement whose
//! *precondition* is unrecoverable reads exactly like one whose precondition was
//! what you assumed.
//!
//! # The stamp
//!
//! [`ErrorModeStamp`] is derived from the [`ErrorConfig`] the generator actually
//! built, never written by hand:
//!
//! ```
//! use ferro_hgvs::conformance::error_mode_stamp::ErrorModeStamp;
//! use ferro_hgvs::error_handling::{ErrorConfig, ErrorOverride, ErrorType};
//!
//! let config = ErrorConfig::lenient().with_override(ErrorType::MissingVersion, ErrorOverride::Reject);
//! let stamp = ErrorModeStamp::of(&config);
//!
//! assert_eq!(stamp.error_mode, "lenient");
//! assert_eq!(stamp.error_mode_overrides["MissingVersion"], "reject");
//! ```
//!
//! Three properties are deliberate, and each is a lesson from `CaptureCounts`
//! being stampable-but-never-stamped:
//!
//! - **It carries the overrides, not only the mode.** A mode is the *base* of a
//!   per-code resolution (see [`ErrorMode`]'s type docs), so `"lenient"` alone
//!   is not a reproducible precondition when the generator overrode two codes.
//!   An artifact that named only the mode would be a more convincing version of
//!   the same gap.
//! - **There is no [`Default`], and no public constructor that invents one.**
//!   `ErrorModeStamp::of(&config)` is the only way to build one, so a stamp
//!   cannot be filled in from what the author believed the mode was — which is
//!   the belief that was wrong in every case above.
//! - **It is [`Deserialize`] as well as [`Serialize`]**, so a consuming test can
//!   assert the precondition it is comparing against rather than assume it.
//!
//! **Name the half the stamp covers.** A stamp describes one `ErrorConfig`, and
//! a generator usually has two stages that could carry one. Both spec generators
//! here parse with the bare [`parse_hgvs`](crate::parse_hgvs), which constructs
//! no `InputPreprocessor` and so applies **no** `ErrorConfig` at all — not
//! lenient's repairs and not strict's refusals (#1632, pinned by
//! `tests/it/issue_1632_parse_entry_applies_no_mode.rs`). "No mode" is a third
//! thing, not a synonym for strict, so a field named `generated_under` on those
//! artifacts would assert a parser-level precondition no row ever had, which is
//! a fresh ambiguity introduced by the field added to remove one. They name the
//! field `normalized_under` for that reason. Stamp what you configured, and let
//! the field name say how far it reaches.
//!
//! The mechanical half — "a generator that normalizes and writes an artifact
//! must stamp this, or name itself in an allowlist" — is
//! `tests/it/generator_completeness.rs`. Read its reach honestly: it is a
//! substring scan, a floor rather than a proof, exactly as the ledger guard
//! beside it is.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::error_handling::{ErrorConfig, ErrorMode};

/// The error-handling precondition of a generated artifact.
///
/// Stamp it beside an artifact's other identity fields (`description`,
/// `spec`, `completeness`) so a reader can tell which mode produced the numbers
/// without re-running the generator. Build it with [`ErrorModeStamp::of`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ErrorModeStamp {
    /// The base mode, as [`ErrorMode`]'s `Display` spells it — `strict`,
    /// `lenient` or `silent`.
    ///
    /// A string rather than the enum so that reading an artifact written by a
    /// future ferro cannot fail on a mode this build does not know; the point of
    /// the stamp is to be legible, and a hard deserialization error on the
    /// provenance field would defeat that.
    pub error_mode: String,
    /// Per-code overrides in force, keyed by [`ErrorType`]'s `Debug` name and
    /// valued by [`ErrorOverride`]'s `Display`.
    ///
    /// Present because warning and correction are resolved **per code**: the
    /// mode is only the base, and `--ignore`/`--reject` move a code off it in
    /// either direction (#1196, #1629). A stamp naming only the mode would
    /// describe a precondition the run did not have.
    ///
    /// Ordered, and omitted from the artifact when empty — the overwhelmingly
    /// common case — so stamping costs an unoverridden generator one line.
    ///
    /// [`ErrorType`]: crate::error_handling::ErrorType
    /// [`ErrorOverride`]: crate::error_handling::ErrorOverride
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub error_mode_overrides: BTreeMap<String, String>,
}

impl ErrorModeStamp {
    /// Derive the stamp from the config a generator actually built.
    ///
    /// The only constructor, on purpose: see the module docs.
    #[must_use]
    pub fn of(config: &ErrorConfig) -> Self {
        Self {
            error_mode: config.mode.to_string(),
            error_mode_overrides: config
                .overrides
                .iter()
                .map(|(error_type, over)| (format!("{error_type:?}"), over.to_string()))
                .collect(),
        }
    }

    /// Whether this stamp records `mode` with no overrides — i.e. the artifact
    /// was produced under a plain preset.
    ///
    /// For a consuming test that wants to state its precondition rather than
    /// assume it: `assert!(doc.normalized_under.is_plain(ErrorMode::Lenient))`
    /// fails both when the mode moved and when an override was introduced, which
    /// a bare string comparison on `error_mode` would miss.
    #[must_use]
    pub fn is_plain(&self, mode: ErrorMode) -> bool {
        self.error_mode == mode.to_string() && self.error_mode_overrides.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_handling::{ErrorOverride, ErrorType};

    #[test]
    fn a_preset_stamps_its_mode_and_no_overrides() {
        for (config, expected) in [
            (ErrorConfig::strict(), ErrorMode::Strict),
            (ErrorConfig::lenient(), ErrorMode::Lenient),
            (ErrorConfig::silent(), ErrorMode::Silent),
        ] {
            let stamp = ErrorModeStamp::of(&config);
            assert_eq!(stamp.error_mode, expected.to_string());
            assert!(stamp.error_mode_overrides.is_empty());
            assert!(stamp.is_plain(expected));
        }
    }

    /// The property the mode alone cannot carry: two configs that a
    /// mode-only stamp would call identical are not the same precondition.
    #[test]
    fn an_override_is_part_of_the_precondition_and_shows_in_the_stamp() {
        let plain = ErrorModeStamp::of(&ErrorConfig::lenient());
        let overridden = ErrorModeStamp::of(
            &ErrorConfig::lenient().with_override(ErrorType::MissingVersion, ErrorOverride::Reject),
        );

        assert_eq!(plain.error_mode, overridden.error_mode);
        assert_ne!(
            plain, overridden,
            "an override must make the stamp differ; otherwise the stamp claims a \
             precondition the run did not have"
        );
        assert!(!overridden.is_plain(ErrorMode::Lenient));
        assert_eq!(
            overridden.error_mode_overrides["MissingVersion"],
            ErrorOverride::Reject.to_string()
        );
    }

    /// `NormalizeConfig::default()` is lenient while `ErrorConfig::default()` is
    /// strict — the trap the stamp exists for. Pinned so that a change to either
    /// default surfaces here rather than as a silently re-based census.
    #[test]
    fn the_stamp_distinguishes_the_two_defaults_that_disagree() {
        assert_eq!(
            ErrorModeStamp::of(&ErrorConfig::default()).error_mode,
            "strict"
        );
        assert_eq!(
            ErrorModeStamp::of(&crate::NormalizeConfig::default().error_config).error_mode,
            "lenient",
            "NormalizeConfig::default() is LENIENT; a generator written as \
             `Normalizer::new(provider)` measures in lenient mode while reading as default"
        );
    }

    #[test]
    fn the_stamp_round_trips_through_json() {
        let stamp = ErrorModeStamp::of(
            &ErrorConfig::silent()
                .with_override(ErrorType::WrongDashCharacter, ErrorOverride::WarnCorrect),
        );
        let json = serde_json::to_string(&stamp).expect("serialize");
        let read: ErrorModeStamp = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(stamp, read);
    }
}
