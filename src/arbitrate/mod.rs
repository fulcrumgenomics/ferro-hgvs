//! HGVS arbitration: explain and adjudicate parse/normalize/projection
//! disagreements between ferro and another tool. See the design spec
//! `specs/2026-07-09-hgvs-arbitration-skill-design.md`.

pub mod bug_report;
pub mod category;
pub mod oracle;
pub mod predicates;
pub mod spec_citations;

pub use category::ArbitrationCategory;

use crate::arbitrate::oracle::{compare_variants, SameVariant};
use crate::arbitrate::spec_citations::{lookup, Operation, SpecCitation};
use crate::error::FerroError;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;
use crate::reference::provider::ReferenceProvider;
use serde::{Deserialize, Serialize};

/// The full arbitration verdict for one input, produced by [`arbitrate`]:
/// ferro's and the other tool's outputs reduced to a same/different verdict,
/// a spec-compliance judgment on genuine disagreements, and the SPDI forms of
/// both outputs for a bug report. See Task 10 (CLI) and Task 11/12
/// (bug-report) for consumers.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Arbitration {
    /// The original query input (as given to the arbitration, not necessarily
    /// identical to `ferro_output`/`other.output`).
    pub input: String,
    /// Same/different/mismatch verdict from the normalizer-independent oracle.
    pub verdict: Verdict,
    /// Which tool (if either) is spec-compliant on a genuine disagreement.
    pub compliance: Compliance,
    /// The overall arbitration category (who is correct), mirroring
    /// `compliance` in the vocabulary the bug-report/CLI consumers expect.
    pub category: ArbitrationCategory,
    /// Ferro's raw output string, when it produced one.
    pub ferro_output: Option<String>,
    /// The other tool's result as supplied by the caller.
    pub other: OtherResult,
    /// Spec passages that govern this verdict/compliance decision.
    pub spec_citations: Vec<SpecCitation>,
    /// Ferro's output rendered as SPDI, for the bug report.
    pub ferro_spdi: Option<String>,
    /// The other tool's output rendered as SPDI, for the bug report.
    pub other_spdi: Option<String>,
}

/// Same/different/mismatch verdict from the normalizer-independent oracle,
/// or a reason the oracle could not even be run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Verdict {
    /// Both outputs denote the same edit to the reference (`SameVariant::Same`).
    /// `compliance`/`category` may still flag a notation-only issue (e.g. one
    /// side spelling a duplication as an insertion) even though the edit
    /// itself is equivalent.
    Equivalent,
    /// Both outputs denote genuinely different edits.
    Different,
    /// The two outputs reduce to SPDI on different accessions/reference
    /// bases, so there is no shared basis to compare on.
    BasisMismatch,
    /// The other tool's output could not be parsed as HGVS (or was absent).
    OtherUnparseable,
    /// Ferro's own output could not be parsed as HGVS.
    FerroParseError,
}

/// Which tool (if either) is spec-compliant for a genuine disagreement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Compliance {
    /// Ferro's output is the spec-compliant one.
    Ferro,
    /// The other tool's output is the spec-compliant one.
    Other,
    /// Neither output is spec-compliant.
    BothWrong,
    /// The disagreement is genuine but no deterministic predicate applies;
    /// needs a human to interpret against the spec.
    NeedsInterpretation,
    /// Compliance does not apply (e.g. the verdict was `Equivalent` or a
    /// parse/basis failure short-circuited the comparison).
    NotApplicable,
}

impl std::fmt::Display for Compliance {
    /// Snake_case, matching this enum's `#[serde(rename_all = "snake_case")]`
    /// wire form — used for user-facing rendering (e.g. `ferro arbitrate
    /// --format text`) so `compliance` prints in the same casing convention
    /// as `category` ([`ArbitrationCategory`]'s `Display`) rather than the
    /// PascalCase `{:?}` Debug form.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Compliance::Ferro => "ferro",
            Compliance::Other => "other",
            Compliance::BothWrong => "both_wrong",
            Compliance::NeedsInterpretation => "needs_interpretation",
            Compliance::NotApplicable => "not_applicable",
        };
        write!(f, "{s}")
    }
}

/// The other tool's result, as supplied by the caller (e.g. a mutalyzer
/// invocation recorded by the CLI or a bug-report harness).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OtherResult {
    /// The tool's name (e.g. `"mutalyzer"`).
    pub tool: String,
    /// The tool's reported status (e.g. `"ok"`, `"error"`).
    pub status: String,
    /// The tool's output string, if it produced one.
    pub output: Option<String>,
}

/// Frame resolver abstraction so `arbitrate` need not name the projector's
/// generic parameter. Implemented for `VariantProjector<P>` in Task 11 wiring.
pub trait FrameResolver {
    /// Resolve `v` into a frame the oracle can apply directly to reference
    /// (e.g. project an intronic transcript variant onto its genomic parent).
    fn resolve(&self, v: &HgvsVariant) -> Result<HgvsVariant, FerroError>;
}

/// Build the "nothing further to compute" `Arbitration` shell shared by every
/// early-exit branch (parse failures, basis mismatch): only `verdict`,
/// `category`, and `compliance` vary; the rest are the input's identity plus
/// empty derived fields.
fn shell(
    input: &str,
    ferro_output: Option<&str>,
    other: OtherResult,
    verdict: Verdict,
    category: ArbitrationCategory,
    compliance: Compliance,
) -> Arbitration {
    Arbitration {
        input: input.to_string(),
        verdict,
        compliance,
        category,
        ferro_output: ferro_output.map(str::to_string),
        other,
        spec_citations: Vec::new(),
        ferro_spdi: None,
        other_spdi: None,
    }
}

/// Classify a [`SameVariant::Same`] verdict: the two outputs edit the
/// reference identically, but one may still violate the dup-vs-ins notation
/// rule (`duplication.md`: a duplication must be spelled as `dup`, not as the
/// equivalent insertion). [`predicates::dup_vs_ins_compliance`] only returns
/// `Some` when both sides reduce to the *same* reference-anchored edit tuple
/// — which is exactly the condition that makes the oracle report `Same` in
/// the first place — so this check belongs here, not on `Different` (a
/// genuine edit disagreement can never satisfy the predicate's identical-edit
/// precondition).
///
/// Returns `(category, compliance, spec_citation_operation)`; the
/// `ArbitrationCategory::Equivalent`/`Compliance::NotApplicable` default is
/// used when the predicate does not apply (no notation issue detected).
fn classify_same(
    ferro: &HgvsVariant,
    other: &HgvsVariant,
    provider: &(impl ReferenceProvider + ?Sized),
) -> (ArbitrationCategory, Compliance, Operation) {
    match crate::arbitrate::predicates::dup_vs_ins_compliance(ferro, other, provider) {
        Some(ArbitrationCategory::FerroCorrect) => (
            ArbitrationCategory::FerroCorrect,
            Compliance::Ferro,
            Operation::DupVsIns,
        ),
        Some(ArbitrationCategory::MutalyzerCorrect) => (
            ArbitrationCategory::MutalyzerCorrect,
            Compliance::Other,
            Operation::DupVsIns,
        ),
        _ => (
            ArbitrationCategory::Equivalent,
            Compliance::NotApplicable,
            Operation::ThreePrimeShift,
        ),
    }
}

/// Orchestrate the full arbitration: parse both outputs, resolve frames if a
/// projector is given, run the normalizer-independent oracle, classify a
/// genuine disagreement via the deterministic compliance predicates, attach
/// governing spec citations, and populate the SPDI forms of both outputs for
/// a bug report.
///
/// `projector` is `None` when only exonic/genomic comparison is needed (the
/// oracle handles those directly); pass a [`FrameResolver`] to also resolve
/// intronic transcript variants onto their genomic parent first.
///
/// # Errors
///
/// Returns [`FerroError`] if frame resolution fails (when a projector is
/// given and either variant is intronic and unplaceable) or if the oracle
/// cannot reduce a variant to SPDI (e.g. an unsupported edit, or the
/// provider lacks the requested reference data). Parse failures on either
/// side are NOT propagated as `Err` — they are reported as
/// [`Verdict::FerroParseError`]/[`Verdict::OtherUnparseable`] in the returned
/// `Arbitration`, since a failure to parse is itself part of the verdict.
pub fn arbitrate(
    input: &str,
    ferro_output: &str,
    other: OtherResult,
    provider: &(impl ReferenceProvider + ?Sized),
    // projector optional: None => exonic/genomic-only comparison
    projector: Option<&dyn FrameResolver>,
) -> Result<Arbitration, FerroError> {
    let vf = match parse_hgvs(ferro_output) {
        Ok(v) => v,
        Err(_) => {
            return Ok(shell(
                input,
                Some(ferro_output),
                other,
                Verdict::FerroParseError,
                ArbitrationCategory::Unknown,
                Compliance::NotApplicable,
            ))
        }
    };
    let other_str = match &other.output {
        Some(s) => s.clone(),
        None => {
            return Ok(shell(
                input,
                Some(ferro_output),
                other,
                Verdict::OtherUnparseable,
                ArbitrationCategory::Unknown,
                Compliance::NotApplicable,
            ))
        }
    };
    let vo = match parse_hgvs(&other_str) {
        Ok(v) => v,
        Err(_) => {
            return Ok(shell(
                input,
                Some(ferro_output),
                other,
                Verdict::OtherUnparseable,
                ArbitrationCategory::Unknown,
                Compliance::NotApplicable,
            ))
        }
    };

    let (vf, vo) = match projector {
        Some(pr) => (pr.resolve(&vf)?, pr.resolve(&vo)?),
        None => (vf, vo),
    };

    let mut out = match compare_variants(&vf, &vo, provider)? {
        SameVariant::Same => {
            let (category, compliance, op) = classify_same(&vf, &vo, provider);
            let mut a = shell(
                input,
                Some(ferro_output),
                other,
                Verdict::Equivalent,
                category,
                compliance,
            );
            a.spec_citations = lookup(op);
            a
        }
        SameVariant::Different => {
            // No deterministic predicate currently covers genuine edit
            // disagreements (predicates.rs only has the identical-edit
            // dup-vs-ins rule, handled above) — flag for human review.
            let mut a = shell(
                input,
                Some(ferro_output),
                other,
                Verdict::Different,
                ArbitrationCategory::Unknown,
                Compliance::NeedsInterpretation,
            );
            a.spec_citations = lookup(Operation::ThreePrimeShift);
            a
        }
        SameVariant::BasisMismatch { .. } => shell(
            input,
            Some(ferro_output),
            other,
            Verdict::BasisMismatch,
            ArbitrationCategory::Unknown,
            Compliance::NotApplicable,
        ),
    };
    out.ferro_spdi = crate::spdi::convert::hgvs_to_spdi(&vf, provider)
        .ok()
        .map(|s| s.to_string());
    out.other_spdi = crate::spdi::convert::hgvs_to_spdi(&vo, provider)
        .ok()
        .map(|s| s.to_string());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::mock::MockProvider;

    /// A 105 bp synthetic genomic contig with an A-free left flank (1-based
    /// 1..=50), a 5 bp poly-A run (1-based 51..=55), and an A-free right
    /// flank (1-based 56..=105). The poly-A run lets a single-base
    /// duplication be spelled at several positions inside the run while
    /// denoting the same variant — mirrors `oracle::verdict_tests`'
    /// `provider_with_a_run`, since `arbitrate` needs the same
    /// shift-equivalence property to exercise `Verdict::Equivalent` and no
    /// `Transcript::for_test` helper exists in this codebase (test data here
    /// is genomic `g.`, not transcript `c.`/`n.`).
    fn provider() -> MockProvider {
        let left: String = "CGT".chars().cycle().take(50).collect();
        let right: String = "GTC".chars().cycle().take(50).collect();
        let seq = format!("{left}AAAAA{right}");
        let mut p = MockProvider::new();
        p.add_genomic_sequence("NC_000001.11", seq);
        p
    }

    #[test]
    fn equivalent_shift_forms_yield_equivalent_verdict() {
        let p = provider();
        let other = OtherResult {
            tool: "mutalyzer".into(),
            status: "ok".into(),
            output: Some("NC_000001.11:g.55dup".into()),
        };
        let a = arbitrate(
            "NC_000001.11:g.51dup",
            "NC_000001.11:g.51dup",
            other,
            &p,
            None,
        )
        .unwrap();
        assert_eq!(a.verdict, Verdict::Equivalent);
        assert_eq!(a.category, ArbitrationCategory::Equivalent);
        assert_eq!(a.compliance, Compliance::NotApplicable);
        assert!(!a.spec_citations.is_empty());
        assert!(a.ferro_spdi.is_some());
        assert!(a.other_spdi.is_some());
    }

    #[test]
    fn other_unparseable_is_flagged_not_crash() {
        let p = provider();
        let other = OtherResult {
            tool: "mutalyzer".into(),
            status: "ok".into(),
            output: Some("not-hgvs".into()),
        };
        let a = arbitrate(
            "NC_000001.11:g.51dup",
            "NC_000001.11:g.51dup",
            other,
            &p,
            None,
        )
        .unwrap();
        assert_eq!(a.verdict, Verdict::OtherUnparseable);
        assert_eq!(a.compliance, Compliance::NotApplicable);
    }

    #[test]
    fn missing_other_output_is_unparseable_not_a_panic() {
        let p = provider();
        let other = OtherResult {
            tool: "mutalyzer".into(),
            status: "error".into(),
            output: None,
        };
        let a = arbitrate(
            "NC_000001.11:g.51dup",
            "NC_000001.11:g.51dup",
            other,
            &p,
            None,
        )
        .unwrap();
        assert_eq!(a.verdict, Verdict::OtherUnparseable);
    }

    #[test]
    fn unparseable_ferro_output_is_flagged() {
        let p = provider();
        let other = OtherResult {
            tool: "mutalyzer".into(),
            status: "ok".into(),
            output: Some("NC_000001.11:g.51dup".into()),
        };
        let a = arbitrate("NC_000001.11:g.51dup", "garbage", other, &p, None).unwrap();
        assert_eq!(a.verdict, Verdict::FerroParseError);
        assert_eq!(a.compliance, Compliance::NotApplicable);
    }

    #[test]
    fn dup_vs_ins_of_the_same_edit_is_equivalent_but_flags_ferro_correct() {
        // 10 bp contig ATG|CG|CG|TAA; g.4_5dup and g.5_6insCG are the same
        // tandem-duplication edit spelled two ways (mirrors
        // `predicates::tests::provider`). The oracle reduces both to an
        // identical reference-anchored edit tuple, so the verdict is
        // `Equivalent` (same edit) — but the dup-vs-ins notation predicate
        // still flags that only ferro spelled it the spec-compliant way.
        let mut p = MockProvider::new();
        p.add_genomic_sequence("NC_000002.11", "ATGCGCGTAA");
        let other = OtherResult {
            tool: "mutalyzer".into(),
            status: "ok".into(),
            output: Some("NC_000002.11:g.5_6insCG".into()),
        };
        let a = arbitrate(
            "NC_000002.11:g.4_5dup",
            "NC_000002.11:g.4_5dup",
            other,
            &p,
            None,
        )
        .unwrap();
        assert_eq!(a.verdict, Verdict::Equivalent);
        assert_eq!(a.category, ArbitrationCategory::FerroCorrect);
        assert_eq!(a.compliance, Compliance::Ferro);
        assert!(!a.spec_citations.is_empty());
    }

    #[test]
    fn different_accessions_are_basis_mismatch() {
        let mut p = provider();
        p.add_genomic_sequence("NC_000002.11", "T".repeat(105));
        let other = OtherResult {
            tool: "mutalyzer".into(),
            status: "ok".into(),
            output: Some("NC_000002.11:g.51dup".into()),
        };
        let a = arbitrate(
            "NC_000001.11:g.51dup",
            "NC_000001.11:g.51dup",
            other,
            &p,
            None,
        )
        .unwrap();
        assert_eq!(a.verdict, Verdict::BasisMismatch);
        assert_eq!(a.compliance, Compliance::NotApplicable);
    }
}
