//! Reference provider trait
//!
//! Defines the interface for accessing reference sequence data.

use crate::error::FerroError;
use crate::hgvs::variant::HgvsVariant;
use crate::reference::transcript::Transcript;

/// Trait for providing reference sequence data
///
/// Implementations might include:
/// - MockProvider for testing
/// - SeqRepoProvider for local sequence databases
/// - NCBIProvider for remote NCBI E-utilities
pub trait ReferenceProvider {
    /// Get a transcript by its accession
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError>;

    /// Get the transcript for the given variant, optionally using the
    /// variant's `genomic_context` parent (NG/NC) to pick a cdot build
    /// that has the transcript with a chromosome populated.
    ///
    /// The default implementation delegates to
    /// [`get_transcript`](Self::get_transcript) using
    /// [`HgvsVariant::accession`] / [`Accession::transcript_accession`],
    /// which preserves existing behavior for inputs that do not carry a
    /// `genomic_context`.
    ///
    /// Providers that know how to resolve an NG/NC-anchored input to a
    /// specific cdot build (e.g. [`crate::reference::multi_fasta::MultiFastaProvider`])
    /// override this to consult the parent accession before falling back to
    /// the bare transcript lookup. See issue #332.
    fn get_transcript_for_variant(&self, variant: &HgvsVariant) -> Result<Transcript, FerroError> {
        let accession = variant
            .accession()
            .ok_or_else(|| FerroError::ReferenceNotFound {
                id: format!("{}", variant),
            })?;
        self.get_transcript(&accession.transcript_accession())
    }

    /// Get a sequence region
    ///
    /// # Arguments
    ///
    /// * `id` - Sequence accession
    /// * `start` - 0-based start position
    /// * `end` - 0-based end position (exclusive)
    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError>;

    /// Check if a transcript exists
    fn has_transcript(&self, id: &str) -> bool {
        self.get_transcript(id).is_ok()
    }

    /// Get genomic sequence for a contig/chromosome
    ///
    /// This is used for normalizing intronic variants which require access
    /// to genomic (not transcript) sequence data.
    ///
    /// # Arguments
    ///
    /// * `contig` - Chromosome/contig accession (e.g., "NC_000001.11", "chr1")
    /// * `start` - 0-based start position
    /// * `end` - 0-based end position (exclusive)
    ///
    /// # Returns
    ///
    /// The genomic sequence, or an error if genomic data is not available.
    /// The default implementation returns an error indicating genomic data
    /// is not available.
    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        Err(FerroError::GenomicReferenceNotAvailable {
            contig: contig.to_string(),
            start,
            end,
        })
    }

    /// Check if this provider has genomic sequence data
    ///
    /// Returns true if `get_genomic_sequence` can return data.
    fn has_genomic_data(&self) -> bool {
        false
    }

    /// Get protein sequence for a protein accession
    ///
    /// This is used for validating protein variants against the reference
    /// protein sequence.
    ///
    /// # Arguments
    ///
    /// * `accession` - Protein accession (e.g., "NP_000079.2", "ENSP00000256509")
    /// * `start` - 0-based start position (amino acid)
    /// * `end` - 0-based end position (exclusive)
    ///
    /// # Returns
    ///
    /// The protein sequence (amino acid letters), or an error if not available.
    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        Err(FerroError::ProteinReferenceNotAvailable {
            accession: accession.to_string(),
            start,
            end,
        })
    }

    /// Check if this provider has protein sequence data
    fn has_protein_data(&self) -> bool {
        false
    }

    /// Get the total length of a stored sequence by accession.
    ///
    /// Used by circular-aware span math (issue #399) to compute the
    /// wraparound length `(L − start + 1) + end` on `m.`/`o.` ranges
    /// where `start > end`. Also consumed by W4004 PositionPastEnd on
    /// the `m.` axis (issue #393).
    ///
    /// The default implementation returns
    /// [`FerroError::ReferenceNotFound`]; providers that store length
    /// metadata override this.
    fn get_seq_length(&self, id: &str) -> Result<u64, FerroError> {
        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }
}

/// Blanket implementation for boxed trait objects
impl ReferenceProvider for Box<dyn ReferenceProvider> {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        (**self).get_transcript(id)
    }

    fn get_transcript_for_variant(&self, variant: &HgvsVariant) -> Result<Transcript, FerroError> {
        (**self).get_transcript_for_variant(variant)
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        (**self).get_sequence(id, start, end)
    }

    fn has_transcript(&self, id: &str) -> bool {
        (**self).has_transcript(id)
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        (**self).get_genomic_sequence(contig, start, end)
    }

    fn has_genomic_data(&self) -> bool {
        (**self).has_genomic_data()
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        (**self).get_protein_sequence(accession, start, end)
    }

    fn has_protein_data(&self) -> bool {
        (**self).has_protein_data()
    }

    fn get_seq_length(&self, id: &str) -> Result<u64, FerroError> {
        (**self).get_seq_length(id)
    }
}
