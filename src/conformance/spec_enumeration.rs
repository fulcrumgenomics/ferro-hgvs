//! Typed classification vocabularies for the spec-enumeration corpus.
//!
//! The exhaustive HGVS spec enumeration (`examples/generate_spec_enumeration.rs`)
//! emits one JSON `Row` per assertion, and the test driver
//! (`tests/it/spec_enumeration_tests.rs`) reads those rows back. Three of a
//! row's classification fields draw from small, closed vocabularies:
//!
//! * [`Status`] — the outcome class (`repaired`, `false-acceptance`, …);
//! * [`Expectation`] — whether the row pins a spec expectation or ferro's
//!   current behaviour;
//! * [`NormativeLevel`] — the RFC 2119 force of the row.
//!
//! Modelling them as enums shared by the generator and the driver turns a typo
//! or a drift between the two — previously caught only at runtime — into a
//! compile error. The generator and driver stay decoupled through the JSON
//! fixture, so the runtime "is every observed status accounted for?" guard from
//! issue #1107 remains valuable: [`Status`] therefore keeps an
//! [`Unknown`](Status::Unknown) fallback so a *new* status the generator starts
//! emitting deserializes and is named by that guard, rather than degrading into
//! an opaque serde "unknown variant" error.

use std::fmt;

use serde::de::Deserializer;
use serde::{Deserialize, Serialize, Serializer};

/// The outcome class recorded for one enumerated assertion.
///
/// The wire form is the kebab-case string the generator has always emitted, so
/// the generated fixture is byte-for-byte unchanged. Any string that is not a
/// known variant deserializes into [`Unknown`](Status::Unknown), preserving the
/// offending name for the #1107 accounted-for diagnostic instead of failing the
/// whole parse.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Status {
    /// Spec-forbidden string that ferro rejected outright (parse error) — the
    /// conformant counterpart of [`FalseAcceptance`](Status::FalseAcceptance).
    CorrectlyRejected,
    /// Spec-forbidden string that ferro parses and renders anyway.
    FalseAcceptance,
    /// Bad string repaired to the canonical form the spec names.
    Repaired,
    /// Spec names a canonical replacement; ferro rejects the bad string instead
    /// of repairing it. Spec-conformant — rejection is always permitted.
    RejectedNotRepaired,
    /// Spec names a canonical replacement; ferro accepts the bad string and
    /// renders something else. A genuine violation.
    RepairDiverges,
    /// Repair that needs real reference bases; not assertable hermetically.
    RequiresReference,
    /// Per-error-mode outcome pinned as ferro policy (the spec says nothing
    /// about ferro's error modes).
    ModeDivergencePinned,
    /// Grammar example that parses into the coordinate axis the spec declares.
    FormAxisOk,
    /// syntax.yaml example that does not parse into its declared axis.
    FormAxisDiverges,
    /// Grammar example whose axis the spec does not state; ferro's parsed axis
    /// is pinned as a baseline.
    FormAxisPinned,
    /// MUST-level output invariant violated by a string ferro emits.
    InvariantViolationMust,
    /// SHOULD-level (advisory) output invariant. Never a hard failure.
    InvariantViolationShould,
    /// Projection that matches the form the spec states — the conformant
    /// counterpart of [`ProjectionDiverges`](Status::ProjectionDiverges).
    Preserved,
    /// Projection whose rendered form differs from the one the spec states.
    ProjectionDiverges,
    /// Projection with no spec-stated expectation, rendered and pinned as a
    /// baseline.
    ProjectionPinned,
    /// Projection with no spec-stated expectation that ferro declined /
    /// reported unavailable, pinned as a baseline.
    ProjectionUnavailablePinned,
    /// Projection with no spec-stated expectation that errored, pinned as a
    /// baseline.
    ProjectionErrorPinned,
    /// A status string the generator emitted that is not (yet) a known variant.
    /// Its name is preserved so the #1107 accounted-for guard can report it.
    Unknown(String),
}

impl Status {
    /// The kebab-case wire string for this status. For
    /// [`Unknown`](Status::Unknown) this is the preserved original name.
    pub fn as_str(&self) -> &str {
        match self {
            Status::CorrectlyRejected => "correctly-rejected",
            Status::FalseAcceptance => "false-acceptance",
            Status::Repaired => "repaired",
            Status::RejectedNotRepaired => "rejected-not-repaired",
            Status::RepairDiverges => "repair-diverges",
            Status::RequiresReference => "requires-reference",
            Status::ModeDivergencePinned => "mode-divergence-pinned",
            Status::FormAxisOk => "form-axis-ok",
            Status::FormAxisDiverges => "form-axis-diverges",
            Status::FormAxisPinned => "form-axis-pinned",
            Status::InvariantViolationMust => "invariant-violation-must",
            Status::InvariantViolationShould => "invariant-violation-should",
            Status::Preserved => "preserved",
            Status::ProjectionDiverges => "projection-diverges",
            Status::ProjectionPinned => "projection-pinned",
            Status::ProjectionUnavailablePinned => "projection-unavailable-pinned",
            Status::ProjectionErrorPinned => "projection-error-pinned",
            Status::Unknown(s) => s,
        }
    }

    /// Parse a wire string into a status, mapping any unrecognised value to
    /// [`Unknown`](Status::Unknown) so the name survives for diagnostics.
    pub fn from_wire(s: &str) -> Self {
        match s {
            "correctly-rejected" => Status::CorrectlyRejected,
            "false-acceptance" => Status::FalseAcceptance,
            "repaired" => Status::Repaired,
            "rejected-not-repaired" => Status::RejectedNotRepaired,
            "repair-diverges" => Status::RepairDiverges,
            "requires-reference" => Status::RequiresReference,
            "mode-divergence-pinned" => Status::ModeDivergencePinned,
            "form-axis-ok" => Status::FormAxisOk,
            "form-axis-diverges" => Status::FormAxisDiverges,
            "form-axis-pinned" => Status::FormAxisPinned,
            "invariant-violation-must" => Status::InvariantViolationMust,
            "invariant-violation-should" => Status::InvariantViolationShould,
            "preserved" => Status::Preserved,
            "projection-diverges" => Status::ProjectionDiverges,
            "projection-pinned" => Status::ProjectionPinned,
            "projection-unavailable-pinned" => Status::ProjectionUnavailablePinned,
            "projection-error-pinned" => Status::ProjectionErrorPinned,
            other => Status::Unknown(other.to_string()),
        }
    }
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl Serialize for Status {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for Status {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        Ok(Status::from_wire(&s))
    }
}

/// Whether a row states a spec expectation or merely pins ferro's current
/// behaviour. A `pinned-baseline` row must never carry MUST-level force.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Expectation {
    /// The spec states the expected value outright.
    SpecMandated,
    /// The spec states nothing; the row records current behaviour instead.
    PinnedBaseline,
}

impl Expectation {
    /// The kebab-case wire string for this expectation.
    pub fn as_str(&self) -> &'static str {
        match self {
            Expectation::SpecMandated => "spec-mandated",
            Expectation::PinnedBaseline => "pinned-baseline",
        }
    }
}

impl fmt::Display for Expectation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// RFC 2119 normative force of a row. Only [`Must`](NormativeLevel::Must) rows
/// may ever hard-fail; [`Na`](NormativeLevel::Na) marks a pinned baseline the
/// spec places no obligation on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum NormativeLevel {
    /// RFC 2119 MUST.
    #[serde(rename = "must")]
    Must,
    /// RFC 2119 SHOULD.
    #[serde(rename = "should")]
    Should,
    /// No normative force (pinned baseline).
    #[serde(rename = "n/a")]
    Na,
}

impl NormativeLevel {
    /// The wire string for this level (`must` / `should` / `n/a`).
    pub fn as_str(&self) -> &'static str {
        match self {
            NormativeLevel::Must => "must",
            NormativeLevel::Should => "should",
            NormativeLevel::Na => "n/a",
        }
    }
}

impl fmt::Display for NormativeLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every known status must round-trip name → variant → name unchanged, and
    /// serialize to exactly that wire string.
    #[test]
    fn status_known_variants_round_trip() {
        let names = [
            "correctly-rejected",
            "false-acceptance",
            "repaired",
            "rejected-not-repaired",
            "repair-diverges",
            "requires-reference",
            "mode-divergence-pinned",
            "form-axis-ok",
            "form-axis-diverges",
            "form-axis-pinned",
            "invariant-violation-must",
            "invariant-violation-should",
            "preserved",
            "projection-diverges",
            "projection-pinned",
            "projection-unavailable-pinned",
            "projection-error-pinned",
        ];
        for name in names {
            let status = Status::from_wire(name);
            assert!(
                !matches!(status, Status::Unknown(_)),
                "{name} should map to a known variant"
            );
            assert_eq!(status.as_str(), name);
            let json = serde_json::to_string(&status).unwrap();
            assert_eq!(json, format!("\"{name}\""));
            let back: Status = serde_json::from_str(&json).unwrap();
            assert_eq!(back, status);
        }
    }

    /// An unknown status deserializes into `Unknown` with its name preserved,
    /// rather than failing the parse — this is what keeps the #1107 diagnostic
    /// able to name the offender.
    #[test]
    fn status_unknown_preserves_name() {
        let status: Status = serde_json::from_str("\"projection panicked\"").unwrap();
        assert_eq!(status, Status::Unknown("projection panicked".to_string()));
        assert_eq!(status.as_str(), "projection panicked");
        // And round-trips back out verbatim.
        assert_eq!(
            serde_json::to_string(&status).unwrap(),
            "\"projection panicked\""
        );
    }

    #[test]
    fn expectation_round_trips() {
        for (variant, wire) in [
            (Expectation::SpecMandated, "\"spec-mandated\""),
            (Expectation::PinnedBaseline, "\"pinned-baseline\""),
        ] {
            assert_eq!(serde_json::to_string(&variant).unwrap(), wire);
            assert_eq!(serde_json::from_str::<Expectation>(wire).unwrap(), variant);
        }
    }

    #[test]
    fn normative_level_round_trips() {
        for (variant, wire) in [
            (NormativeLevel::Must, "\"must\""),
            (NormativeLevel::Should, "\"should\""),
            (NormativeLevel::Na, "\"n/a\""),
        ] {
            assert_eq!(serde_json::to_string(&variant).unwrap(), wire);
            assert_eq!(
                serde_json::from_str::<NormativeLevel>(wire).unwrap(),
                variant
            );
        }
    }

    /// An unknown expectation / normative level fails loudly (naming the value),
    /// since these small closed vocabularies have no decoupled runtime guard.
    #[test]
    fn strict_enums_reject_unknown() {
        assert!(serde_json::from_str::<Expectation>("\"made-up\"").is_err());
        assert!(serde_json::from_str::<NormativeLevel>("\"maybe\"").is_err());
    }
}
