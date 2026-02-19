//! ClinVar API client implementation.
//!
//! This module provides a client for accessing ClinVar data. Currently implements
//! a mock client for testing; real API access can be added with an HTTP client.

use std::collections::HashMap;

use super::types::ClinVarRecord;

/// ClinVar API client.
///
/// This client provides access to ClinVar variant data. The mock implementation
/// stores records in memory for testing purposes.
///
/// # Example
///
/// ```
/// use ferro_hgvs::clinvar::{ClinVarClient, ClinVarRecord, ClinicalSignificance};
///
/// let mut client = ClinVarClient::mock();
///
/// // Add a test record
/// client.add_record(ClinVarRecord {
///     variation_id: "12345".to_string(),
///     hgvs: "NM_000088.3:c.10A>G".to_string(),
///     significance: ClinicalSignificance::Pathogenic,
///     review_status: "criteria provided, single submitter".to_string(),
///     ..Default::default()
/// });
///
/// // Query by variation ID
/// let record = client.get_by_variation_id("12345");
/// assert!(record.is_some());
/// ```
#[derive(Debug, Clone)]
pub struct ClinVarClient {
    /// Records indexed by variation ID.
    records_by_id: HashMap<String, ClinVarRecord>,
    /// Records indexed by HGVS expression.
    records_by_hgvs: HashMap<String, String>, // HGVS -> variation_id
    /// Records indexed by rsID.
    records_by_rsid: HashMap<String, String>, // rsID -> variation_id
    /// Records indexed by accession.
    records_by_accession: HashMap<String, String>, // accession -> variation_id
}

impl Default for ClinVarClient {
    fn default() -> Self {
        Self::new()
    }
}

impl ClinVarClient {
    /// Create a new empty client.
    pub fn new() -> Self {
        Self {
            records_by_id: HashMap::new(),
            records_by_hgvs: HashMap::new(),
            records_by_rsid: HashMap::new(),
            records_by_accession: HashMap::new(),
        }
    }

    /// Create a mock client for testing.
    ///
    /// This is equivalent to `new()` but makes the intent clearer.
    pub fn mock() -> Self {
        Self::new()
    }

    /// Add a record to the client.
    ///
    /// This indexes the record by variation ID, HGVS, rsID (if present),
    /// and accession (if present).
    pub fn add_record(&mut self, record: ClinVarRecord) {
        let id = record.variation_id.clone();

        // Index by HGVS
        self.records_by_hgvs.insert(record.hgvs.clone(), id.clone());

        // Index by HGVS aliases
        for alias in &record.hgvs_aliases {
            self.records_by_hgvs.insert(alias.clone(), id.clone());
        }

        // Index by rsID if present
        if let Some(ref rsid) = record.rsid {
            self.records_by_rsid.insert(rsid.clone(), id.clone());
        }

        // Index by accession if present
        if let Some(ref accession) = record.accession {
            self.records_by_accession
                .insert(accession.clone(), id.clone());
        }

        // Store the record
        self.records_by_id.insert(id, record);
    }

    /// Get a record by ClinVar variation ID.
    pub fn get_by_variation_id(&self, id: &str) -> Option<&ClinVarRecord> {
        self.records_by_id.get(id)
    }

    /// Get a record by HGVS expression.
    ///
    /// This searches both the primary HGVS and any aliases.
    pub fn get_by_hgvs(&self, hgvs: &str) -> Option<&ClinVarRecord> {
        self.records_by_hgvs
            .get(hgvs)
            .and_then(|id| self.records_by_id.get(id))
    }

    /// Get a record by rsID.
    pub fn get_by_rsid(&self, rsid: &str) -> Option<&ClinVarRecord> {
        // Normalize rsID format
        let normalized = if rsid.starts_with("rs") {
            rsid.to_string()
        } else {
            format!("rs{}", rsid)
        };

        self.records_by_rsid
            .get(&normalized)
            .or_else(|| self.records_by_rsid.get(rsid))
            .and_then(|id| self.records_by_id.get(id))
    }

    /// Get a record by ClinVar accession (e.g., "VCV000012345").
    pub fn get_by_accession(&self, accession: &str) -> Option<&ClinVarRecord> {
        self.records_by_accession
            .get(accession)
            .and_then(|id| self.records_by_id.get(id))
    }

    /// Search for records by gene symbol.
    ///
    /// Returns all records associated with the given gene.
    pub fn search_by_gene(&self, gene: &str) -> Vec<&ClinVarRecord> {
        let gene_lower = gene.to_lowercase();
        self.records_by_id
            .values()
            .filter(|r| {
                r.gene
                    .as_ref()
                    .map(|g| g.to_lowercase() == gene_lower)
                    .unwrap_or(false)
            })
            .collect()
    }

    /// Search for pathogenic variants in a gene.
    ///
    /// Returns records classified as Pathogenic or Likely Pathogenic.
    pub fn search_pathogenic_in_gene(&self, gene: &str) -> Vec<&ClinVarRecord> {
        self.search_by_gene(gene)
            .into_iter()
            .filter(|r| r.significance.is_pathogenic())
            .collect()
    }

    /// Get all records in the client.
    pub fn all_records(&self) -> impl Iterator<Item = &ClinVarRecord> {
        self.records_by_id.values()
    }

    /// Get the number of records in the client.
    pub fn len(&self) -> usize {
        self.records_by_id.len()
    }

    /// Check if the client has no records.
    pub fn is_empty(&self) -> bool {
        self.records_by_id.is_empty()
    }

    /// Remove a record by variation ID.
    ///
    /// Returns the removed record if it existed.
    pub fn remove(&mut self, variation_id: &str) -> Option<ClinVarRecord> {
        if let Some(record) = self.records_by_id.remove(variation_id) {
            // Clean up indexes
            self.records_by_hgvs.remove(&record.hgvs);
            for alias in &record.hgvs_aliases {
                self.records_by_hgvs.remove(alias);
            }
            if let Some(ref rsid) = record.rsid {
                self.records_by_rsid.remove(rsid);
            }
            if let Some(ref accession) = record.accession {
                self.records_by_accession.remove(accession);
            }
            Some(record)
        } else {
            None
        }
    }

    /// Clear all records from the client.
    pub fn clear(&mut self) {
        self.records_by_id.clear();
        self.records_by_hgvs.clear();
        self.records_by_rsid.clear();
        self.records_by_accession.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clinvar::{ClinicalSignificance, ConditionInfo, ReviewStatus};

    fn create_test_record() -> ClinVarRecord {
        ClinVarRecord {
            variation_id: "12345".to_string(),
            accession: Some("VCV000012345".to_string()),
            hgvs: "NM_000088.3:c.10A>G".to_string(),
            hgvs_aliases: vec!["NC_000017.11:g.48275363A>G".to_string()],
            gene: Some("COL1A1".to_string()),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "criteria provided, multiple submitters, no conflicts".to_string(),
            parsed_review_status: Some(ReviewStatus::MultipleSubmitters),
            conditions: vec![ConditionInfo::new("Osteogenesis imperfecta")],
            rsid: Some("rs12345".to_string()),
            ..Default::default()
        }
    }

    #[test]
    fn test_new_client() {
        let client = ClinVarClient::new();
        assert!(client.is_empty());
        assert_eq!(client.len(), 0);
    }

    #[test]
    fn test_mock_client() {
        let client = ClinVarClient::mock();
        assert!(client.is_empty());
    }

    #[test]
    fn test_add_and_get_by_id() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        let record = client.get_by_variation_id("12345");
        assert!(record.is_some());
        assert_eq!(record.unwrap().variation_id, "12345");
    }

    #[test]
    fn test_get_by_hgvs() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        // Primary HGVS
        let record = client.get_by_hgvs("NM_000088.3:c.10A>G");
        assert!(record.is_some());

        // Alias HGVS
        let record = client.get_by_hgvs("NC_000017.11:g.48275363A>G");
        assert!(record.is_some());

        // Non-existent
        let record = client.get_by_hgvs("NM_999999.1:c.1A>G");
        assert!(record.is_none());
    }

    #[test]
    fn test_get_by_rsid() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        // With rs prefix
        let record = client.get_by_rsid("rs12345");
        assert!(record.is_some());

        // Without rs prefix (should still work)
        let record = client.get_by_rsid("12345");
        assert!(record.is_some());
    }

    #[test]
    fn test_get_by_accession() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        let record = client.get_by_accession("VCV000012345");
        assert!(record.is_some());
        assert_eq!(record.unwrap().variation_id, "12345");
    }

    #[test]
    fn test_search_by_gene() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        // Add another record for a different gene
        let mut other_record = create_test_record();
        other_record.variation_id = "67890".to_string();
        other_record.gene = Some("BRCA1".to_string());
        other_record.hgvs = "NM_007294.4:c.100A>G".to_string();
        other_record.hgvs_aliases = vec![];
        other_record.rsid = None;
        other_record.accession = None;
        client.add_record(other_record);

        // Search for COL1A1
        let results = client.search_by_gene("COL1A1");
        assert_eq!(results.len(), 1);

        // Case insensitive
        let results = client.search_by_gene("col1a1");
        assert_eq!(results.len(), 1);

        // Search for BRCA1
        let results = client.search_by_gene("BRCA1");
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn test_search_pathogenic_in_gene() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        // Add a benign variant
        let mut benign = create_test_record();
        benign.variation_id = "11111".to_string();
        benign.significance = ClinicalSignificance::Benign;
        benign.hgvs = "NM_000088.3:c.20A>G".to_string();
        benign.hgvs_aliases = vec![];
        benign.rsid = None;
        benign.accession = None;
        client.add_record(benign);

        let pathogenic = client.search_pathogenic_in_gene("COL1A1");
        assert_eq!(pathogenic.len(), 1);
        assert!(pathogenic[0].significance.is_pathogenic());
    }

    #[test]
    fn test_all_records() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        let all: Vec<_> = client.all_records().collect();
        assert_eq!(all.len(), 1);
    }

    #[test]
    fn test_remove() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        assert_eq!(client.len(), 1);

        let removed = client.remove("12345");
        assert!(removed.is_some());
        assert_eq!(client.len(), 0);

        // All indexes should be cleaned
        assert!(client.get_by_hgvs("NM_000088.3:c.10A>G").is_none());
        assert!(client.get_by_rsid("rs12345").is_none());
        assert!(client.get_by_accession("VCV000012345").is_none());
    }

    #[test]
    fn test_clear() {
        let mut client = ClinVarClient::mock();
        client.add_record(create_test_record());

        client.clear();
        assert!(client.is_empty());
    }

    // Fixture-based tests

    /// Create a client populated with realistic ClinVar data
    fn create_fixture_client() -> ClinVarClient {
        let mut client = ClinVarClient::mock();

        // TP53 variants (tumor suppressor)
        client.add_record(ClinVarRecord {
            variation_id: "12375".to_string(),
            accession: Some("VCV000012375".to_string()),
            hgvs: "NM_000546.6:c.215C>G".to_string(),
            hgvs_aliases: vec!["NC_000017.11:g.7674220C>G".to_string()],
            gene: Some("TP53".to_string()),
            significance: ClinicalSignificance::LikelyBenign,
            review_status: "criteria provided, multiple submitters, no conflicts".to_string(),
            parsed_review_status: Some(ReviewStatus::MultipleSubmitters),
            conditions: vec![ConditionInfo::new("Li-Fraumeni syndrome")],
            rsid: Some("rs1042522".to_string()),
            ..Default::default()
        });

        client.add_record(ClinVarRecord {
            variation_id: "12356".to_string(),
            accession: Some("VCV000012356".to_string()),
            hgvs: "NM_000546.6:c.743G>A".to_string(),
            hgvs_aliases: vec!["NC_000017.11:g.7673803G>A".to_string()],
            gene: Some("TP53".to_string()),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "reviewed by expert panel".to_string(),
            parsed_review_status: Some(ReviewStatus::ExpertPanel),
            conditions: vec![
                ConditionInfo::new("Li-Fraumeni syndrome"),
                ConditionInfo::new("Hereditary cancer-predisposing syndrome"),
            ],
            rsid: Some("rs28934576".to_string()),
            ..Default::default()
        });

        // BRCA1 variants (breast cancer)
        client.add_record(ClinVarRecord {
            variation_id: "17661".to_string(),
            accession: Some("VCV000017661".to_string()),
            hgvs: "NM_007294.4:c.5266dupC".to_string(),
            hgvs_aliases: vec![
                "NC_000017.11:g.43057062dupC".to_string(),
                "NP_009225.1:p.Gln1756ProfsTer74".to_string(),
            ],
            gene: Some("BRCA1".to_string()),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "reviewed by expert panel".to_string(),
            parsed_review_status: Some(ReviewStatus::ExpertPanel),
            conditions: vec![ConditionInfo::new(
                "Hereditary breast and ovarian cancer syndrome",
            )],
            rsid: Some("rs80357906".to_string()),
            ..Default::default()
        });

        client.add_record(ClinVarRecord {
            variation_id: "17668".to_string(),
            accession: Some("VCV000017668".to_string()),
            hgvs: "NM_007294.4:c.68_69del".to_string(),
            hgvs_aliases: vec![],
            gene: Some("BRCA1".to_string()),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "reviewed by expert panel".to_string(),
            parsed_review_status: Some(ReviewStatus::ExpertPanel),
            conditions: vec![ConditionInfo::new(
                "Hereditary breast and ovarian cancer syndrome",
            )],
            rsid: Some("rs80357914".to_string()),
            ..Default::default()
        });

        // BRCA2 variants
        client.add_record(ClinVarRecord {
            variation_id: "51557".to_string(),
            accession: Some("VCV000051557".to_string()),
            hgvs: "NM_000059.4:c.5946del".to_string(),
            hgvs_aliases: vec!["NC_000013.11:g.32340300del".to_string()],
            gene: Some("BRCA2".to_string()),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "reviewed by expert panel".to_string(),
            parsed_review_status: Some(ReviewStatus::ExpertPanel),
            conditions: vec![ConditionInfo::new(
                "Hereditary breast and ovarian cancer syndrome",
            )],
            rsid: Some("rs80359550".to_string()),
            ..Default::default()
        });

        // CFTR variant (cystic fibrosis)
        client.add_record(ClinVarRecord {
            variation_id: "7105".to_string(),
            accession: Some("VCV000007105".to_string()),
            hgvs: "NM_000492.4:c.1521_1523del".to_string(),
            hgvs_aliases: vec![
                "NC_000007.14:g.117559590_117559592del".to_string(),
                "NP_000483.3:p.Phe508del".to_string(),
            ],
            gene: Some("CFTR".to_string()),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "reviewed by expert panel".to_string(),
            parsed_review_status: Some(ReviewStatus::ExpertPanel),
            conditions: vec![ConditionInfo::new("Cystic fibrosis")],
            rsid: Some("rs113993960".to_string()),
            ..Default::default()
        });

        // MLH1 variant (Lynch syndrome)
        client.add_record(ClinVarRecord {
            variation_id: "89986".to_string(),
            accession: Some("VCV000089986".to_string()),
            hgvs: "NM_000249.4:c.1852_1854del".to_string(),
            hgvs_aliases: vec![],
            gene: Some("MLH1".to_string()),
            significance: ClinicalSignificance::Pathogenic,
            review_status: "reviewed by expert panel".to_string(),
            parsed_review_status: Some(ReviewStatus::ExpertPanel),
            conditions: vec![ConditionInfo::new("Lynch syndrome")],
            rsid: Some("rs63750447".to_string()),
            ..Default::default()
        });

        // VUS example
        client.add_record(ClinVarRecord {
            variation_id: "12377".to_string(),
            accession: Some("VCV000012377".to_string()),
            hgvs: "NM_000546.6:c.467G>A".to_string(),
            hgvs_aliases: vec![],
            gene: Some("TP53".to_string()),
            significance: ClinicalSignificance::UncertainSignificance,
            review_status: "criteria provided, single submitter".to_string(),
            parsed_review_status: Some(ReviewStatus::SingleSubmitter),
            conditions: vec![ConditionInfo::new("not provided")],
            rsid: None,
            ..Default::default()
        });

        // Benign variant
        client.add_record(ClinVarRecord {
            variation_id: "12400".to_string(),
            accession: Some("VCV000012400".to_string()),
            hgvs: "NM_000546.6:c.639A>G".to_string(),
            hgvs_aliases: vec![],
            gene: Some("TP53".to_string()),
            significance: ClinicalSignificance::Benign,
            review_status: "criteria provided, multiple submitters, no conflicts".to_string(),
            parsed_review_status: Some(ReviewStatus::MultipleSubmitters),
            conditions: vec![ConditionInfo::new("not specified")],
            rsid: Some("rs1800370".to_string()),
            ..Default::default()
        });

        // Conflicting interpretations
        client.add_record(ClinVarRecord {
            variation_id: "99999".to_string(),
            accession: Some("VCV000099999".to_string()),
            hgvs: "NM_000059.4:c.316+5G>A".to_string(),
            hgvs_aliases: vec![],
            gene: Some("BRCA2".to_string()),
            significance: ClinicalSignificance::Conflicting,
            review_status: "criteria provided, conflicting classifications".to_string(),
            parsed_review_status: Some(ReviewStatus::ConflictingInterpretations),
            conditions: vec![],
            rsid: None,
            ..Default::default()
        });

        client
    }

    #[test]
    fn test_fixture_load_multiple_records() {
        let client = create_fixture_client();
        assert_eq!(client.len(), 10);
    }

    #[test]
    fn test_fixture_query_tp53_variants() {
        let client = create_fixture_client();
        let tp53_variants = client.search_by_gene("TP53");
        assert_eq!(tp53_variants.len(), 4); // 2 pathogenic + 1 VUS + 1 benign
    }

    #[test]
    fn test_fixture_query_brca1_pathogenic() {
        let client = create_fixture_client();
        let pathogenic = client.search_pathogenic_in_gene("BRCA1");
        assert_eq!(pathogenic.len(), 2);
        assert!(pathogenic.iter().all(|r| r.significance.is_pathogenic()));
    }

    #[test]
    fn test_fixture_query_brca2_pathogenic() {
        let client = create_fixture_client();
        let pathogenic = client.search_pathogenic_in_gene("BRCA2");
        // Only 1 because the conflicting one is not pathogenic
        assert_eq!(pathogenic.len(), 1);
    }

    #[test]
    fn test_fixture_founder_mutations_by_rsid() {
        let client = create_fixture_client();

        // BRCA1 5382insC (Ashkenazi Jewish founder)
        let record = client.get_by_rsid("rs80357906");
        assert!(record.is_some());
        assert_eq!(record.unwrap().gene.as_deref(), Some("BRCA1"));

        // BRCA1 185delAG (Ashkenazi Jewish founder)
        let record = client.get_by_rsid("rs80357914");
        assert!(record.is_some());
        assert_eq!(record.unwrap().hgvs, "NM_007294.4:c.68_69del");

        // BRCA2 6174delT (Ashkenazi Jewish founder)
        let record = client.get_by_rsid("rs80359550");
        assert!(record.is_some());
        assert_eq!(record.unwrap().gene.as_deref(), Some("BRCA2"));
    }

    #[test]
    fn test_fixture_cftr_f508del() {
        let client = create_fixture_client();

        // Most common CF mutation
        let record = client.get_by_hgvs("NM_000492.4:c.1521_1523del");
        assert!(record.is_some());
        let rec = record.unwrap();
        assert_eq!(rec.gene.as_deref(), Some("CFTR"));
        assert!(rec.significance.is_pathogenic());

        // Can also query by protein notation alias
        let record = client.get_by_hgvs("NP_000483.3:p.Phe508del");
        assert!(record.is_some());
    }

    #[test]
    fn test_fixture_expert_reviewed_variants() {
        let client = create_fixture_client();

        let expert_reviewed: Vec<_> = client
            .all_records()
            .filter(|r| r.parsed_review_status == Some(ReviewStatus::ExpertPanel))
            .collect();

        // 6 variants reviewed by expert panel
        assert_eq!(expert_reviewed.len(), 6);
    }

    #[test]
    fn test_fixture_vus_variants() {
        let client = create_fixture_client();

        let vus: Vec<_> = client
            .all_records()
            .filter(|r| r.significance == ClinicalSignificance::UncertainSignificance)
            .collect();

        assert_eq!(vus.len(), 1);
        assert_eq!(vus[0].variation_id, "12377");
    }

    #[test]
    fn test_fixture_conflicting_variants() {
        let client = create_fixture_client();

        let conflicting: Vec<_> = client
            .all_records()
            .filter(|r| r.significance == ClinicalSignificance::Conflicting)
            .collect();

        assert_eq!(conflicting.len(), 1);
        assert_eq!(conflicting[0].hgvs, "NM_000059.4:c.316+5G>A");
    }

    #[test]
    fn test_fixture_query_by_accession_vcv() {
        let client = create_fixture_client();

        // Query using VCV accession
        let record = client.get_by_accession("VCV000007105");
        assert!(record.is_some());
        assert_eq!(record.unwrap().gene.as_deref(), Some("CFTR"));
    }

    #[test]
    fn test_fixture_genomic_hgvs_aliases() {
        let client = create_fixture_client();

        // Query using genomic coordinates instead of transcript
        let record = client.get_by_hgvs("NC_000017.11:g.7673803G>A");
        assert!(record.is_some());
        let rec = record.unwrap();
        assert_eq!(rec.variation_id, "12356");
        assert_eq!(rec.gene.as_deref(), Some("TP53"));
    }

    #[test]
    fn test_fixture_significance_distribution() {
        let client = create_fixture_client();

        let mut pathogenic = 0;
        let mut benign = 0;
        let mut vus = 0;
        let mut conflicting = 0;

        for record in client.all_records() {
            match record.significance {
                ClinicalSignificance::Pathogenic | ClinicalSignificance::LikelyPathogenic => {
                    pathogenic += 1
                }
                ClinicalSignificance::Benign | ClinicalSignificance::LikelyBenign => benign += 1,
                ClinicalSignificance::UncertainSignificance => vus += 1,
                ClinicalSignificance::Conflicting => conflicting += 1,
                _ => {}
            }
        }

        assert_eq!(pathogenic, 6); // TP53 hotspot, 2 BRCA1, BRCA2, CFTR, MLH1
        assert_eq!(benign, 2); // TP53 polymorphism + benign
        assert_eq!(vus, 1);
        assert_eq!(conflicting, 1);
    }

    #[test]
    fn test_fixture_cancer_gene_panel() {
        let client = create_fixture_client();

        // Simulate querying a cancer gene panel
        let cancer_genes = vec!["TP53", "BRCA1", "BRCA2", "MLH1"];
        let mut panel_pathogenic: Vec<&ClinVarRecord> = vec![];

        for gene in &cancer_genes {
            panel_pathogenic.extend(client.search_pathogenic_in_gene(gene));
        }

        // Should find: TP53 R248Q, 2 BRCA1, 1 BRCA2, 1 MLH1
        assert_eq!(panel_pathogenic.len(), 5);
    }

    #[test]
    fn test_fixture_multiple_conditions() {
        let client = create_fixture_client();

        // TP53 R248Q is associated with multiple conditions
        let record = client.get_by_variation_id("12356");
        assert!(record.is_some());
        let rec = record.unwrap();
        assert_eq!(rec.conditions.len(), 2);
    }

    #[test]
    fn test_fixture_review_status_hierarchy() {
        let client = create_fixture_client();

        // Expert panel > multiple submitters > single submitter > conflicting
        let expert_count = client
            .all_records()
            .filter(|r| r.parsed_review_status == Some(ReviewStatus::ExpertPanel))
            .count();

        let multiple_count = client
            .all_records()
            .filter(|r| r.parsed_review_status == Some(ReviewStatus::MultipleSubmitters))
            .count();

        let single_count = client
            .all_records()
            .filter(|r| r.parsed_review_status == Some(ReviewStatus::SingleSubmitter))
            .count();

        let conflicting_count = client
            .all_records()
            .filter(|r| r.parsed_review_status == Some(ReviewStatus::ConflictingInterpretations))
            .count();

        assert_eq!(expert_count, 6);
        assert_eq!(multiple_count, 2);
        assert_eq!(single_count, 1);
        assert_eq!(conflicting_count, 1);
    }
}
