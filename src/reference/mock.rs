//! Mock reference provider for testing

use crate::error::FerroError;
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::Transcript;
use std::collections::HashMap;
use std::path::Path;
use std::sync::OnceLock;

/// Mock reference provider that loads transcripts from JSON
#[derive(Clone)]
pub struct MockProvider {
    transcripts: HashMap<String, Transcript>,
    proteins: HashMap<String, String>,
    /// Genomic sequences keyed by contig name
    genomic_sequences: HashMap<String, String>,
}

impl MockProvider {
    /// Create an empty mock provider
    pub fn new() -> Self {
        Self {
            transcripts: HashMap::new(),
            proteins: HashMap::new(),
            genomic_sequences: HashMap::new(),
        }
    }

    /// Load transcripts from a JSON file
    pub fn from_json(path: &Path) -> Result<Self, FerroError> {
        let content = std::fs::read_to_string(path)?;
        let transcripts: Vec<Transcript> = serde_json::from_str(&content)?;

        let map: HashMap<String, Transcript> = transcripts
            .into_iter()
            .map(|tx| (tx.id.clone(), tx))
            .collect();

        Ok(Self {
            transcripts: map,
            proteins: HashMap::new(),
            genomic_sequences: HashMap::new(),
        })
    }

    /// Add a transcript to the provider
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.insert(transcript.id.clone(), transcript);
    }

    /// Add a protein sequence to the provider
    pub fn add_protein(&mut self, accession: impl Into<String>, sequence: impl Into<String>) {
        self.proteins.insert(accession.into(), sequence.into());
    }

    /// Add a genomic sequence for a contig/chromosome
    pub fn add_genomic_sequence(&mut self, contig: impl Into<String>, sequence: impl Into<String>) {
        self.genomic_sequences
            .insert(contig.into(), sequence.into());
    }

    /// Create a provider with some test transcripts
    pub fn with_test_data() -> Self {
        use crate::reference::transcript::{Exon, ManeStatus, Strand};

        let mut provider = Self::new();

        // Add a typical coding transcript (MANE Select)
        provider.add_transcript(Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCCCAAGGTGCTGCCCCAGATGCTGCCAGTGCTGCTGCTGCTGCTGCTGCTGCTGCTG".to_string(),
            cds_start: Some(1),
            cds_end: Some(60),
            exons: vec![
                Exon::new(1, 1, 20),
                Exon::new(2, 21, 40),
                Exon::new(3, 41, 60),
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        });

        // Add a transcript with UTRs
        provider.add_transcript(Transcript {
            id: "NM_001234.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: "AAAAATGCCCAAGGGGGGGGGGGGGGGGGGGGGGGGGTAAAAAA".to_string(),
            cds_start: Some(5),
            cds_end: Some(38),
            exons: vec![
                Exon::new(1, 1, 15),
                Exon::new(2, 16, 30),
                Exon::new(3, 31, 44),
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        });

        // Add a minus strand transcript
        provider.add_transcript(Transcript {
            id: "NM_999999.1".to_string(),
            gene_symbol: Some("MINUS".to_string()),
            strand: Strand::Minus,
            sequence: "ATGCATGCATGCATGCATGCATGCATGCATGC".to_string(),
            cds_start: Some(1),
            cds_end: Some(30),
            exons: vec![Exon::new(1, 1, 32)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        });

        // Add a transcript for testing 3' normalization of duplications
        // Sequence has "GAA" at c.8-c.9 to test that c.8dup -> c.9dup (2 A's)
        // and "GTTT" at c.21-c.24 to test that c.22dup -> c.24dup (3 T's)
        // Positions: 1234567890123456789012345678901234567890
        // Sequence:  ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        //                  ^^           ^^^
        //                  c.8-9 (AA)   c.22-24 (TTT)
        provider.add_transcript(Transcript {
            id: "NM_888888.1".to_string(),
            gene_symbol: Some("DUPTEST".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT".to_string(),
            cds_start: Some(1),
            cds_end: Some(39),
            exons: vec![Exon::new(1, 1, 39)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        });

        // Add test protein sequences for validation testing
        // NP_000079.2 is used for BRAF V600E-style testing
        // We need Val at position 600, Arg at 97 for frameshift tests
        // Building a sequence with known positions:
        // - Position 1: M (Met)
        // - Position 97: R (Arg) for frameshift test
        // - Position 600: V (Val) for V600E test
        // - Position 23-25: K, A, E for range tests
        let mut seq = String::new();
        // Positions 1-96 (96 chars): just padding with Ala
        seq.push('M'); // Position 1
        for _ in 2..23 {
            seq.push('A'); // Positions 2-22
        }
        seq.push('K'); // Position 23 - for Lys23
        seq.push('A'); // Position 24 - for middle of range
        seq.push('E'); // Position 25 - for Glu25
        for _ in 26..97 {
            seq.push('A'); // Positions 26-96
        }
        seq.push('R'); // Position 97 - for Arg97 frameshift
        for _ in 98..600 {
            seq.push('A'); // Positions 98-599
        }
        seq.push('V'); // Position 600 - for Val600
        for _ in 601..=700 {
            seq.push('A'); // Positions 601-700
        }
        provider.add_protein("NP_000079.2", seq);

        // A simpler test protein for easier position testing
        // Position 1=M, 2=V, 3=L, 4=S, 5=P, etc.
        provider.add_protein("NP_TEST.1", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFL");

        provider
    }

    /// Get the number of transcripts
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Check if provider is empty
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Get all transcript IDs
    pub fn transcript_ids(&self) -> Vec<&str> {
        self.transcripts.keys().map(|s| s.as_str()).collect()
    }
}

impl Default for MockProvider {
    fn default() -> Self {
        Self::new()
    }
}

impl ReferenceProvider for MockProvider {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        // Handle versioned and unversioned lookups
        if let Some(tx) = self.transcripts.get(id) {
            return Ok(tx.clone());
        }

        // Try without version
        let base_id = id.split('.').next().unwrap_or(id);
        for (key, tx) in &self.transcripts {
            if key.starts_with(base_id) {
                return Ok(tx.clone());
            }
        }

        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        let transcript = self.get_transcript(id)?;

        transcript
            .get_sequence(start, end)
            .map(|s| s.to_string())
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: format!("Position {}-{} out of range for {}", start, end, id),
            })
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Handle versioned and unversioned lookups
        let protein_seq = if let Some(seq) = self.proteins.get(accession) {
            seq.clone()
        } else {
            // Try without version
            let base_id = accession.split('.').next().unwrap_or(accession);
            let mut found = None;
            for (key, seq) in &self.proteins {
                if key.starts_with(base_id) {
                    found = Some(seq.clone());
                    break;
                }
            }
            found.ok_or_else(|| FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start,
                end,
            })?
        };

        let start = start as usize;
        let end = end as usize;

        if start >= protein_seq.len() || end > protein_seq.len() || start > end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}-{} out of range for {} (length {})",
                    start,
                    end,
                    accession,
                    protein_seq.len()
                ),
            });
        }

        Ok(protein_seq[start..end].to_string())
    }

    fn has_protein_data(&self) -> bool {
        !self.proteins.is_empty()
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let genomic_seq = self.genomic_sequences.get(contig).ok_or_else(|| {
            FerroError::GenomicReferenceNotAvailable {
                contig: contig.to_string(),
                start,
                end,
            }
        })?;

        let start = start as usize;
        let end = end as usize;

        if start >= genomic_seq.len() || end > genomic_seq.len() || start > end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}-{} out of range for {} (length {})",
                    start,
                    end,
                    contig,
                    genomic_seq.len()
                ),
            });
        }

        Ok(genomic_seq[start..end].to_string())
    }

    fn has_genomic_data(&self) -> bool {
        !self.genomic_sequences.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mock_provider_with_test_data() {
        let provider = MockProvider::with_test_data();
        assert!(!provider.is_empty());
        assert!(provider.len() >= 2);
    }

    #[test]
    fn test_get_transcript() {
        let provider = MockProvider::with_test_data();
        let tx = provider.get_transcript("NM_000088.3").unwrap();
        assert_eq!(tx.gene_symbol, Some("COL1A1".to_string()));
    }

    #[test]
    fn test_get_transcript_not_found() {
        let provider = MockProvider::with_test_data();
        let result = provider.get_transcript("NM_NONEXISTENT.1");
        assert!(result.is_err());
    }

    #[test]
    fn test_get_sequence() {
        let provider = MockProvider::with_test_data();
        let seq = provider.get_sequence("NM_000088.3", 0, 3).unwrap();
        assert_eq!(seq, "ATG");
    }

    #[test]
    fn test_has_transcript() {
        let provider = MockProvider::with_test_data();
        assert!(provider.has_transcript("NM_000088.3"));
        assert!(!provider.has_transcript("NONEXISTENT"));
    }
}
