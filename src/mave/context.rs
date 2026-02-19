//! MAVE-HGVS parsing context.
//!
//! Provides the target sequence context needed to interpret MAVE-HGVS short forms.

use serde::{Deserialize, Serialize};

/// Context for parsing MAVE-HGVS short-form notation.
///
/// MaveDB score sets define a target sequence that provides the context for
/// interpreting HGVS variants without explicit accessions.
///
/// # Example
///
/// ```
/// use ferro_hgvs::mave::MaveContext;
///
/// let context = MaveContext::new()
///     .with_protein_accession("NP_000509.1")
///     .with_coding_accession("NM_000518.5")
///     .with_gene_symbol("HBB");
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct MaveContext {
    /// Protein sequence accession (for p. variants).
    pub protein_accession: Option<String>,

    /// Coding sequence accession (for c. variants).
    pub coding_accession: Option<String>,

    /// Non-coding transcript accession (for n. variants).
    pub noncoding_accession: Option<String>,

    /// Genomic sequence accession (for g. variants).
    pub genomic_accession: Option<String>,

    /// Gene symbol (informational).
    pub gene_symbol: Option<String>,

    /// MaveDB score set URN.
    pub score_set_urn: Option<String>,

    /// Target sequence name (from MaveDB).
    pub target_name: Option<String>,
}

impl MaveContext {
    /// Create a new empty context.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the protein sequence accession for p. variants.
    pub fn with_protein_accession(mut self, accession: impl Into<String>) -> Self {
        self.protein_accession = Some(accession.into());
        self
    }

    /// Set the coding sequence accession for c. variants.
    pub fn with_coding_accession(mut self, accession: impl Into<String>) -> Self {
        self.coding_accession = Some(accession.into());
        self
    }

    /// Set the non-coding transcript accession for n. variants.
    pub fn with_noncoding_accession(mut self, accession: impl Into<String>) -> Self {
        self.noncoding_accession = Some(accession.into());
        self
    }

    /// Set the genomic sequence accession for g. variants.
    pub fn with_genomic_accession(mut self, accession: impl Into<String>) -> Self {
        self.genomic_accession = Some(accession.into());
        self
    }

    /// Set the gene symbol (informational).
    pub fn with_gene_symbol(mut self, symbol: impl Into<String>) -> Self {
        self.gene_symbol = Some(symbol.into());
        self
    }

    /// Set the MaveDB score set URN.
    pub fn with_score_set_urn(mut self, urn: impl Into<String>) -> Self {
        self.score_set_urn = Some(urn.into());
        self
    }

    /// Set the target sequence name.
    pub fn with_target_name(mut self, name: impl Into<String>) -> Self {
        self.target_name = Some(name.into());
        self
    }

    /// Get the appropriate accession for a coordinate type.
    ///
    /// # Arguments
    ///
    /// * `coord_type` - The coordinate type character ('p', 'c', 'n', 'g', etc.)
    ///
    /// # Returns
    ///
    /// The accession string if available for the given coordinate type.
    pub fn accession_for_coordinate_type(&self, coord_type: char) -> Option<&str> {
        match coord_type {
            'p' => self.protein_accession.as_deref(),
            'c' => self.coding_accession.as_deref(),
            'n' => self.noncoding_accession.as_deref(),
            'g' => self.genomic_accession.as_deref(),
            'r' => self.coding_accession.as_deref(), // RNA uses same as coding
            'm' | 'o' => self.genomic_accession.as_deref(), // mitochondrial uses genomic
            _ => None,
        }
    }

    /// Check if this context has any accessions defined.
    pub fn has_accessions(&self) -> bool {
        self.protein_accession.is_some()
            || self.coding_accession.is_some()
            || self.noncoding_accession.is_some()
            || self.genomic_accession.is_some()
    }

    /// Check if this context can handle a specific coordinate type.
    pub fn supports_coordinate_type(&self, coord_type: char) -> bool {
        self.accession_for_coordinate_type(coord_type).is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_context() {
        let ctx = MaveContext::new();
        assert!(ctx.protein_accession.is_none());
        assert!(ctx.coding_accession.is_none());
        assert!(!ctx.has_accessions());
    }

    #[test]
    fn test_with_protein_accession() {
        let ctx = MaveContext::new().with_protein_accession("NP_000509.1");
        assert_eq!(ctx.protein_accession, Some("NP_000509.1".to_string()));
        assert!(ctx.has_accessions());
    }

    #[test]
    fn test_with_coding_accession() {
        let ctx = MaveContext::new().with_coding_accession("NM_000518.5");
        assert_eq!(ctx.coding_accession, Some("NM_000518.5".to_string()));
        assert!(ctx.has_accessions());
    }

    #[test]
    fn test_full_context() {
        let ctx = MaveContext::new()
            .with_protein_accession("NP_000509.1")
            .with_coding_accession("NM_000518.5")
            .with_genomic_accession("NC_000011.10")
            .with_gene_symbol("HBB")
            .with_score_set_urn("urn:mavedb:00000001-a-1");

        assert_eq!(ctx.protein_accession, Some("NP_000509.1".to_string()));
        assert_eq!(ctx.coding_accession, Some("NM_000518.5".to_string()));
        assert_eq!(ctx.genomic_accession, Some("NC_000011.10".to_string()));
        assert_eq!(ctx.gene_symbol, Some("HBB".to_string()));
        assert_eq!(
            ctx.score_set_urn,
            Some("urn:mavedb:00000001-a-1".to_string())
        );
    }

    #[test]
    fn test_accession_for_coordinate_type() {
        let ctx = MaveContext::new()
            .with_protein_accession("NP_000509.1")
            .with_coding_accession("NM_000518.5")
            .with_genomic_accession("NC_000011.10");

        assert_eq!(ctx.accession_for_coordinate_type('p'), Some("NP_000509.1"));
        assert_eq!(ctx.accession_for_coordinate_type('c'), Some("NM_000518.5"));
        assert_eq!(ctx.accession_for_coordinate_type('g'), Some("NC_000011.10"));
        assert_eq!(ctx.accession_for_coordinate_type('n'), None);
        assert_eq!(ctx.accession_for_coordinate_type('r'), Some("NM_000518.5"));
        // RNA uses coding
    }

    #[test]
    fn test_supports_coordinate_type() {
        let ctx = MaveContext::new()
            .with_protein_accession("NP_000509.1")
            .with_coding_accession("NM_000518.5");

        assert!(ctx.supports_coordinate_type('p'));
        assert!(ctx.supports_coordinate_type('c'));
        assert!(!ctx.supports_coordinate_type('g'));
        assert!(!ctx.supports_coordinate_type('n'));
    }

    #[test]
    fn test_serde() {
        let ctx = MaveContext::new()
            .with_protein_accession("NP_000509.1")
            .with_gene_symbol("HBB");

        let json = serde_json::to_string(&ctx).unwrap();
        let ctx2: MaveContext = serde_json::from_str(&json).unwrap();
        assert_eq!(ctx, ctx2);
    }
}
