//! ClinVar variant information client.
//!
//! This module provides types and utilities for working with ClinVar data.
//! It includes data structures for ClinVar records, a mock client for testing,
//! and integration points for real API access.
//!
//! # Examples
//!
//! ## Using Mock Client for Testing
//!
//! ```
//! use ferro_hgvs::clinvar::{ClinVarClient, ClinVarRecord, ClinicalSignificance};
//!
//! let mut client = ClinVarClient::mock();
//!
//! // Add a test record
//! client.add_record(ClinVarRecord {
//!     variation_id: "12345".to_string(),
//!     hgvs: "NM_000088.3:c.10A>G".to_string(),
//!     significance: ClinicalSignificance::Pathogenic,
//!     review_status: "criteria_provided, single submitter".to_string(),
//!     ..Default::default()
//! });
//!
//! // Query by HGVS
//! let record = client.get_by_hgvs("NM_000088.3:c.10A>G");
//! assert!(record.is_some());
//! ```
//!
//! # API Integration
//!
//! To use with the real ClinVar API (requires HTTP client):
//!
//! ```ignore
//! // Future implementation with reqwest
//! let client = ClinVarClient::new("https://api.ncbi.nlm.nih.gov/clinvar/");
//! let record = client.get_by_variation_id("12345").await?;
//! ```
//!
//! # References
//!
//! - [ClinVar API](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)
//! - [ClinVar Data Dictionary](https://www.ncbi.nlm.nih.gov/clinvar/docs/help/)

mod client;
mod types;

pub use client::ClinVarClient;
pub use types::{ClinVarRecord, ClinicalSignificance, ConditionInfo, ReviewStatus, SubmitterInfo};
