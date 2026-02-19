//! ClinVar data types.

use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// Clinical significance classification from ClinVar.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum ClinicalSignificance {
    /// Benign - variant does not cause disease
    Benign,
    /// Likely benign - variant probably does not cause disease
    LikelyBenign,
    /// Uncertain significance - insufficient evidence
    UncertainSignificance,
    /// Likely pathogenic - variant probably causes disease
    LikelyPathogenic,
    /// Pathogenic - variant causes disease
    Pathogenic,
    /// Conflicting interpretations from different submitters
    Conflicting,
    /// Drug response
    DrugResponse,
    /// Association - variant is associated with phenotype
    Association,
    /// Risk factor
    RiskFactor,
    /// Protective factor
    Protective,
    /// Affects (affects gene function)
    Affects,
    /// Not provided
    #[default]
    NotProvided,
    /// Other
    Other,
}

impl ClinicalSignificance {
    /// Convert to ClinVar string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Benign => "Benign",
            Self::LikelyBenign => "Likely benign",
            Self::UncertainSignificance => "Uncertain significance",
            Self::LikelyPathogenic => "Likely pathogenic",
            Self::Pathogenic => "Pathogenic",
            Self::Conflicting => "Conflicting interpretations of pathogenicity",
            Self::DrugResponse => "drug response",
            Self::Association => "association",
            Self::RiskFactor => "risk factor",
            Self::Protective => "Protective",
            Self::Affects => "Affects",
            Self::NotProvided => "not provided",
            Self::Other => "other",
        }
    }

    /// Check if this is a pathogenic classification.
    pub fn is_pathogenic(&self) -> bool {
        matches!(self, Self::Pathogenic | Self::LikelyPathogenic)
    }

    /// Check if this is a benign classification.
    pub fn is_benign(&self) -> bool {
        matches!(self, Self::Benign | Self::LikelyBenign)
    }

    /// Check if this is uncertain or conflicting.
    pub fn is_uncertain(&self) -> bool {
        matches!(self, Self::UncertainSignificance | Self::Conflicting)
    }
}

impl std::fmt::Display for ClinicalSignificance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl FromStr for ClinicalSignificance {
    type Err = std::convert::Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s_lower = s.to_lowercase();
        Ok(match s_lower.as_str() {
            "benign" => Self::Benign,
            "likely benign" | "likely_benign" => Self::LikelyBenign,
            "uncertain significance" | "uncertain_significance" | "vus" => {
                Self::UncertainSignificance
            }
            "likely pathogenic" | "likely_pathogenic" => Self::LikelyPathogenic,
            "pathogenic" => Self::Pathogenic,
            "conflicting interpretations of pathogenicity"
            | "conflicting_interpretations_of_pathogenicity"
            | "conflicting" => Self::Conflicting,
            "drug response" | "drug_response" => Self::DrugResponse,
            "association" => Self::Association,
            "risk factor" | "risk_factor" => Self::RiskFactor,
            "protective" => Self::Protective,
            "affects" => Self::Affects,
            "not provided" | "not_provided" => Self::NotProvided,
            _ => Self::Other,
        })
    }
}

/// ClinVar review status (star rating).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum ReviewStatus {
    /// No assertion criteria provided (0 stars)
    #[default]
    NoAssertion,
    /// No assertion provided (0 stars)
    NoAssertionProvided,
    /// Criteria provided, conflicting interpretations (1 star)
    ConflictingInterpretations,
    /// Criteria provided, single submitter (1 star)
    SingleSubmitter,
    /// Criteria provided, multiple submitters (2 stars)
    MultipleSubmitters,
    /// Reviewed by expert panel (3 stars)
    ExpertPanel,
    /// Practice guideline (4 stars)
    PracticeGuideline,
}

impl ReviewStatus {
    /// Get the star rating (0-4).
    pub fn stars(&self) -> u8 {
        match self {
            Self::NoAssertion | Self::NoAssertionProvided => 0,
            Self::ConflictingInterpretations | Self::SingleSubmitter => 1,
            Self::MultipleSubmitters => 2,
            Self::ExpertPanel => 3,
            Self::PracticeGuideline => 4,
        }
    }

    /// Convert to ClinVar string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::NoAssertion => "no assertion criteria provided",
            Self::NoAssertionProvided => "no assertion provided",
            Self::ConflictingInterpretations => "criteria provided, conflicting interpretations",
            Self::SingleSubmitter => "criteria provided, single submitter",
            Self::MultipleSubmitters => "criteria provided, multiple submitters, no conflicts",
            Self::ExpertPanel => "reviewed by expert panel",
            Self::PracticeGuideline => "practice guideline",
        }
    }
}

impl std::fmt::Display for ReviewStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl FromStr for ReviewStatus {
    type Err = std::convert::Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s_lower = s.to_lowercase();
        Ok(if s_lower.contains("practice guideline") {
            Self::PracticeGuideline
        } else if s_lower.contains("expert panel") {
            Self::ExpertPanel
        } else if s_lower.contains("multiple submitters") {
            Self::MultipleSubmitters
        } else if s_lower.contains("single submitter") {
            Self::SingleSubmitter
        } else if s_lower.contains("conflicting") {
            Self::ConflictingInterpretations
        } else if s_lower.contains("no assertion provided") {
            Self::NoAssertionProvided
        } else {
            Self::NoAssertion
        })
    }
}

/// Information about a disease/condition.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub struct ConditionInfo {
    /// Condition name.
    pub name: String,
    /// MedGen concept ID (e.g., "C0001080").
    pub medgen_id: Option<String>,
    /// OMIM ID.
    pub omim_id: Option<String>,
    /// Orphanet ID.
    pub orphanet_id: Option<String>,
}

impl ConditionInfo {
    /// Create a new condition with just a name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }

    /// Set the MedGen ID.
    pub fn with_medgen(mut self, id: impl Into<String>) -> Self {
        self.medgen_id = Some(id.into());
        self
    }

    /// Set the OMIM ID.
    pub fn with_omim(mut self, id: impl Into<String>) -> Self {
        self.omim_id = Some(id.into());
        self
    }
}

/// Information about a submitter organization.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub struct SubmitterInfo {
    /// Organization name.
    pub name: String,
    /// Organization ID in ClinVar.
    pub org_id: Option<String>,
    /// Submission date.
    pub date: Option<String>,
}

impl SubmitterInfo {
    /// Create a new submitter info.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }
}

/// A ClinVar variant record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
pub struct ClinVarRecord {
    /// ClinVar Variation ID.
    pub variation_id: String,
    /// ClinVar Accession (e.g., "VCV000012345").
    pub accession: Option<String>,
    /// HGVS expression (typically the canonical expression).
    pub hgvs: String,
    /// Additional HGVS expressions for this variant.
    pub hgvs_aliases: Vec<String>,
    /// Gene symbol.
    pub gene: Option<String>,
    /// Clinical significance.
    pub significance: ClinicalSignificance,
    /// Review status (raw string from ClinVar).
    pub review_status: String,
    /// Parsed review status.
    pub parsed_review_status: Option<ReviewStatus>,
    /// Associated conditions.
    pub conditions: Vec<ConditionInfo>,
    /// Submitter information.
    pub submitters: Vec<SubmitterInfo>,
    /// rsID if available.
    pub rsid: Option<String>,
    /// Last evaluated date.
    pub last_evaluated: Option<String>,
    /// Last updated date in ClinVar.
    pub last_updated: Option<String>,
}

impl ClinVarRecord {
    /// Create a new ClinVar record.
    pub fn new(variation_id: impl Into<String>, hgvs: impl Into<String>) -> Self {
        Self {
            variation_id: variation_id.into(),
            hgvs: hgvs.into(),
            ..Default::default()
        }
    }

    /// Get the review status stars (0-4).
    pub fn stars(&self) -> u8 {
        self.parsed_review_status
            .map(|r| r.stars())
            .unwrap_or_else(|| {
                self.review_status
                    .parse::<ReviewStatus>()
                    .unwrap_or_default()
                    .stars()
            })
    }

    /// Check if this is a high-confidence pathogenic variant (2+ stars).
    pub fn is_high_confidence_pathogenic(&self) -> bool {
        self.significance.is_pathogenic() && self.stars() >= 2
    }

    /// Check if this is a high-confidence benign variant (2+ stars).
    pub fn is_high_confidence_benign(&self) -> bool {
        self.significance.is_benign() && self.stars() >= 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clinical_significance_from_str() {
        assert_eq!(
            "Pathogenic".parse::<ClinicalSignificance>().unwrap(),
            ClinicalSignificance::Pathogenic
        );
        assert_eq!(
            "likely benign".parse::<ClinicalSignificance>().unwrap(),
            ClinicalSignificance::LikelyBenign
        );
        assert_eq!(
            "VUS".parse::<ClinicalSignificance>().unwrap(),
            ClinicalSignificance::UncertainSignificance
        );
        assert_eq!(
            "unknown".parse::<ClinicalSignificance>().unwrap(),
            ClinicalSignificance::Other
        );
    }

    #[test]
    fn test_clinical_significance_is_pathogenic() {
        assert!(ClinicalSignificance::Pathogenic.is_pathogenic());
        assert!(ClinicalSignificance::LikelyPathogenic.is_pathogenic());
        assert!(!ClinicalSignificance::Benign.is_pathogenic());
        assert!(!ClinicalSignificance::UncertainSignificance.is_pathogenic());
    }

    #[test]
    fn test_clinical_significance_is_benign() {
        assert!(ClinicalSignificance::Benign.is_benign());
        assert!(ClinicalSignificance::LikelyBenign.is_benign());
        assert!(!ClinicalSignificance::Pathogenic.is_benign());
    }

    #[test]
    fn test_clinical_significance_is_uncertain() {
        assert!(ClinicalSignificance::UncertainSignificance.is_uncertain());
        assert!(ClinicalSignificance::Conflicting.is_uncertain());
        assert!(!ClinicalSignificance::Pathogenic.is_uncertain());
    }

    #[test]
    fn test_review_status_from_str() {
        assert_eq!(
            "criteria provided, single submitter"
                .parse::<ReviewStatus>()
                .unwrap(),
            ReviewStatus::SingleSubmitter
        );
        assert_eq!(
            "reviewed by expert panel".parse::<ReviewStatus>().unwrap(),
            ReviewStatus::ExpertPanel
        );
        assert_eq!(
            "practice guideline".parse::<ReviewStatus>().unwrap(),
            ReviewStatus::PracticeGuideline
        );
    }

    #[test]
    fn test_review_status_stars() {
        assert_eq!(ReviewStatus::NoAssertion.stars(), 0);
        assert_eq!(ReviewStatus::SingleSubmitter.stars(), 1);
        assert_eq!(ReviewStatus::MultipleSubmitters.stars(), 2);
        assert_eq!(ReviewStatus::ExpertPanel.stars(), 3);
        assert_eq!(ReviewStatus::PracticeGuideline.stars(), 4);
    }

    #[test]
    fn test_clinvar_record_new() {
        let record = ClinVarRecord::new("12345", "NM_000088.3:c.10A>G");
        assert_eq!(record.variation_id, "12345");
        assert_eq!(record.hgvs, "NM_000088.3:c.10A>G");
    }

    #[test]
    fn test_clinvar_record_stars() {
        let record = ClinVarRecord {
            variation_id: "12345".to_string(),
            hgvs: "NM_000088.3:c.10A>G".to_string(),
            review_status: "criteria provided, single submitter".to_string(),
            ..Default::default()
        };
        assert_eq!(record.stars(), 1);
    }

    #[test]
    fn test_clinvar_record_high_confidence_pathogenic() {
        let record = ClinVarRecord {
            variation_id: "12345".to_string(),
            hgvs: "NM_000088.3:c.10A>G".to_string(),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "criteria provided, multiple submitters".to_string(),
            ..Default::default()
        };
        assert!(record.is_high_confidence_pathogenic());

        let low_confidence = ClinVarRecord {
            review_status: "no assertion criteria provided".to_string(),
            ..record.clone()
        };
        assert!(!low_confidence.is_high_confidence_pathogenic());
    }

    #[test]
    fn test_condition_info_builder() {
        let condition = ConditionInfo::new("Hereditary cancer syndrome")
            .with_medgen("C1234567")
            .with_omim("123456");

        assert_eq!(condition.name, "Hereditary cancer syndrome");
        assert_eq!(condition.medgen_id, Some("C1234567".to_string()));
        assert_eq!(condition.omim_id, Some("123456".to_string()));
    }
}
