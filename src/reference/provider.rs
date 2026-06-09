//! Reference provider trait
//!
//! Defines the interface for accessing reference sequence data.

use std::sync::Arc;

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
    /// Get a transcript by its accession.
    ///
    /// Returns an [`Arc`] so that providers which materialize the full
    /// transcript sequence (e.g.
    /// [`MultiFastaProvider`](crate::reference::multi_fasta::MultiFastaProvider))
    /// can memoize the result and hand out cheap pointer clones on repeated
    /// lookups of the same transcript. Callers that only read fields are
    /// unaffected (`Arc<Transcript>` derefs to `Transcript`); callers that
    /// need an owned `Transcript` can clone via `(*tx).clone()`.
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError>;

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
    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
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

    /// Whether `id`'s bases are available at the *exact* requested versioned
    /// accession — i.e. not via a silent version-strip substitution to a
    /// sibling version.
    ///
    /// Protein prediction reads a transcript's CDS bases directly and
    /// translates them. If the reference only carries a *different* version of
    /// the transcript, those bases may pair with a different CDS/exon
    /// structure, so the prediction would be wrong yet attributed to the
    /// requested version (issue #505). The protein path consults this before
    /// the CDS read and declines to predict when it returns `false`.
    ///
    /// The default is `true`: providers that do not track multiple versions of
    /// a transcript (e.g. test doubles) never trigger the gate, so existing
    /// behavior is preserved. Version-aware providers (e.g.
    /// `MultiFastaProvider`) override this to mean "the exact version is
    /// directly present, or can be reconstructed at the exact version" — which
    /// deliberately still permits exact-version cdot-genome synthesis and only
    /// rejects cross-version substitution.
    fn has_transcript_version_exact(&self, _id: &str) -> bool {
        true
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

    /// Return the length (in amino acids) of the named protein.
    ///
    /// # Arguments
    ///
    /// * `accession` - Protein accession (e.g., `"NP_000079.2"`)
    ///
    /// # Returns
    ///
    /// The protein length in amino acids, or `0` if the accession cannot be
    /// resolved by this provider — matching the historical probe's semantics
    /// (see below), so a missing protein behaves identically to before. (The
    /// `Result` is retained for providers that surface genuine I/O errors.)
    ///
    /// The default implementation discovers the length by probing
    /// [`get_protein_sequence`](Self::get_protein_sequence) with
    /// `get_protein_sequence(accession, 0, n)` and binary-searching for
    /// the largest accepted `n`, which equals the protein length. This
    /// preserves the exact semantics of the historical probe: a protein
    /// of length zero (or one the provider cannot resolve) yields a
    /// length of `0`. Providers that store length metadata override this
    /// to return the length directly without cloning the sequence.
    fn get_protein_length(&self, accession: &str) -> Result<u64, FerroError> {
        const SEED: u64 = 64 * 1024;
        const CAP: u64 = 1 << 30;

        let probe_ok = |n: u64| self.get_protein_sequence(accession, 0, n).is_ok();

        let (mut lo, mut hi) = if probe_ok(SEED) {
            // Common case: the seed succeeded. Grow only if the protein
            // is even longer (vanishingly rare).
            let mut hi = SEED;
            let mut lo = SEED;
            while probe_ok(hi) {
                lo = hi;
                if hi >= CAP {
                    break;
                }
                hi = hi.saturating_mul(2);
            }
            (lo, hi)
        } else {
            // Seed failed, so the exact length is in `[0, SEED)`.
            (0, SEED)
        };
        while lo + 1 < hi {
            let mid = lo + (hi - lo) / 2;
            if probe_ok(mid) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    }

    /// Check if this provider has protein sequence data
    fn has_protein_data(&self) -> bool {
        false
    }

    /// Return the total number of bases in the named reference sequence.
    ///
    /// # Arguments
    ///
    /// * `id` - Sequence accession or contig name (e.g., `"NC_012920.1"`, `"chrM"`)
    ///
    /// # Returns
    ///
    /// The length in bases, or [`FerroError::ReferenceNotFound`] if the
    /// accession is unknown to this provider.
    ///
    /// The default implementation always returns
    /// [`FerroError::ReferenceNotFound`]; providers that store length
    /// metadata (e.g. from a FASTA index) override this method.
    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }
}

/// Blanket implementation for boxed trait objects
impl ReferenceProvider for Box<dyn ReferenceProvider> {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        (**self).get_transcript(id)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        (**self).get_transcript_for_variant(variant)
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        (**self).get_sequence(id, start, end)
    }

    fn has_transcript(&self, id: &str) -> bool {
        (**self).has_transcript(id)
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        (**self).has_transcript_version_exact(id)
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

    fn get_protein_length(&self, accession: &str) -> Result<u64, FerroError> {
        (**self).get_protein_length(accession)
    }

    fn has_protein_data(&self) -> bool {
        (**self).has_protein_data()
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        (**self).get_sequence_length(id)
    }
}
