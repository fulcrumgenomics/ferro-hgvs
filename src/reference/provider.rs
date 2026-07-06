//! Reference provider trait
//!
//! Defines the interface for accessing reference sequence data.

use std::sync::Arc;

use crate::error::FerroError;
use crate::hgvs::variant::{Accession, HgvsVariant};
use crate::reference::transcript::Transcript;
use crate::reference::Strand;

/// Chromosomal placement of a genomic *parent* reference — an `NG_` RefSeqGene
/// or an `LRG_` — on its `NC_` chromosome contig, as a single contiguous affine
/// span.
///
/// cdot aligns transcripts only to chromosomes (`NM_`→`NC_`), so a c./n./r.
/// input parented on an `NG_`/`LRG_` resolves to chromosome coordinates that do
/// not match the parent accession's own frame. This placement lets those
/// chromosome coordinates be re-expressed in the parent's frame (#480).
///
/// All coordinates are 1-based. `parent_start` is the parent-frame coordinate of
/// the low-strand span endpoint — the base at `nc_start` on the plus strand, or
/// the base at `nc_end` on the minus strand (see [`GenomicPlacement::nc_to_parent`]);
/// `nc_start`/`nc_end` are the inclusive chromosome span
/// (`nc_start <= nc_end`); `strand` is the parent's orientation relative to the
/// chromosome (`Minus` ⇒ the parent sequence is the reverse complement of the
/// chromosome region, so coordinates count down as chromosome coordinates count
/// up, and edits must be reverse-complemented into the parent frame).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicPlacement {
    /// The chromosome (`NC_`) accession the parent is placed on.
    pub nc: Accession,
    /// Parent-frame (1-based) coordinate of the low-strand span endpoint: the
    /// base at `nc_start` on the plus strand, or the base at `nc_end` on the
    /// minus strand (the parent's own first base). See
    /// [`GenomicPlacement::nc_to_parent`].
    pub parent_start: u64,
    /// Inclusive 1-based chromosome start of the placed span.
    pub nc_start: u64,
    /// Inclusive 1-based chromosome end of the placed span (`>= nc_start`).
    pub nc_end: u64,
    /// Parent orientation relative to the chromosome.
    pub strand: Strand,
}

/// Select the placement matching a resolved genome build from a parent's
/// per-build placement list (#653).
///
/// - `Some(build)` → the placement whose `nc` is on that build, or `None`
///   (decline) when absent — never substitute another build's placement.
/// - `None` (no preference) → prefer the GRCh38 placement (cdot's primary build
///   / mutalyzer policy), else the first available.
pub(crate) fn select_placement_for_build(
    placements: &[GenomicPlacement],
    build: Option<&str>,
    infer_build: impl Fn(&Accession) -> Option<&'static str>,
) -> Option<GenomicPlacement> {
    let build_of = |p: &GenomicPlacement| infer_build(&p.nc);
    match build {
        Some(b) => placements.iter().find(|p| build_of(p) == Some(b)).cloned(),
        None => placements
            .iter()
            .find(|p| build_of(p) == Some("GRCh38"))
            .or_else(|| placements.first())
            .cloned(),
    }
}

impl GenomicPlacement {
    /// Map a 1-based chromosome (`NC_`) coordinate into the parent's own frame.
    ///
    /// Returns `None` when `nc_pos` lies outside the placed span, or when the
    /// parent's strand is unknown — the caller should then decline to re-anchor
    /// rather than emit an out-of-frame coordinate.
    pub fn nc_to_parent(&self, nc_pos: u64) -> Option<u64> {
        if nc_pos < self.nc_start || nc_pos > self.nc_end {
            return None;
        }
        match self.strand {
            Strand::Plus => Some(self.parent_start + (nc_pos - self.nc_start)),
            Strand::Minus => Some(self.parent_start + (self.nc_end - nc_pos)),
            // An unknown-orientation parent has no defined parent frame, so
            // decline rather than guess `Plus`. This matches how the rest of the
            // codebase treats `Strand::Unknown` (transcript projection, annotation
            // validation), and the placement parsers only ever emit Plus/Minus.
            Strand::Unknown => None,
        }
    }

    /// Map a 1-based parent-frame coordinate back onto its chromosome (`NC_`)
    /// coordinate — the inverse of [`Self::nc_to_parent`] (#480). Used to take an
    /// `NG_`/`LRG_`-relative genomic *input* position into the chromosome frame so
    /// cdot can enumerate the transcripts overlapping it.
    ///
    /// Returns `None` when `parent_pos` lies outside the placed span (parent
    /// coordinates `parent_start ..= parent_start + (nc_end - nc_start)`).
    pub fn parent_to_nc(&self, parent_pos: u64) -> Option<u64> {
        let parent_end = self.parent_start + (self.nc_end - self.nc_start);
        if parent_pos < self.parent_start || parent_pos > parent_end {
            return None;
        }
        let offset = parent_pos - self.parent_start;
        match self.strand {
            Strand::Plus | Strand::Unknown => Some(self.nc_start + offset),
            Strand::Minus => Some(self.nc_end - offset),
        }
    }
}

/// Trait for providing reference sequence data
///
/// Implementations might include:
/// - MockProvider for testing
/// - SeqRepoProvider for local sequence databases
/// - NCBIProvider for remote NCBI E-utilities
///
/// **Maintainers:** any method added here with a default implementation must
/// also be forwarded in the `Arc<T>` and `Box<T>` blanket impls below. Rust
/// gives no compile-time guarantee of this, so an un-forwarded method silently
/// falls through to its default when a provider is wrapped (the bug fixed for
/// `resolve_legacy_gene_selector` and `get_transcript_for_accession`).
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
        self.get_transcript_for_accession(accession)
    }

    /// Get the transcript for the given accession, using its `genomic_context`
    /// parent (NG/NC) to pick a cdot build that has the transcript with a
    /// chromosome populated.
    ///
    /// This is the build-aware resolution that backs
    /// [`get_transcript_for_variant`](Self::get_transcript_for_variant), exposed
    /// directly so callers that already hold an [`Accession`] (e.g. the intronic
    /// normalization path) can resolve without cloning a whole variant just to
    /// satisfy the by-variant signature.
    ///
    /// The default implementation delegates to
    /// [`get_transcript`](Self::get_transcript), preserving existing behavior
    /// for inputs that do not carry a `genomic_context`. Build-aware providers
    /// (e.g. [`crate::reference::multi_fasta::MultiFastaProvider`]) override this
    /// to consult the parent accession first. See issue #332.
    fn get_transcript_for_accession(
        &self,
        accession: &Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        self.get_transcript(&accession.transcript_accession())
    }

    /// Return the chromosomal placement of a genomic *parent* reference
    /// (`NG_`/`LRG_`) so transcript coordinates that cdot resolves on the
    /// chromosome can be re-expressed in the parent's own frame (#480).
    ///
    /// The default returns `None` (no placement known); providers that ingest
    /// RefSeqGene / LRG placement data override this. A `None` result means the
    /// projector keeps the chromosome coordinates as-is.
    fn genomic_placement(&self, _parent: &Accession) -> Option<GenomicPlacement> {
        None
    }

    /// Return the chromosomal placement of a genomic *parent* reference for a
    /// specific genome build (#653). `build` is `"GRCh37"`/`"GRCh38"` (or `None`
    /// for "no preference", which prefers GRCh38). When an explicit build is
    /// requested but no placement on that build exists, returns `None` (decline)
    /// rather than mis-anchor onto a different build's placement.
    ///
    /// The default ignores the build hint and delegates to [`genomic_placement`]
    /// (build-agnostic); providers that store per-build placements override this.
    fn genomic_placement_on_build(
        &self,
        parent: &Accession,
        _build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        self.genomic_placement(parent)
    }

    /// Infer the genome build (`"GRCh37"`/`"GRCh38"`) for a chromosome
    /// accession. The default is the bundled hardcoded heuristic; providers
    /// that carry an assembly-report-derived table override this to consult it
    /// first (#716).
    fn infer_genome_build(&self, accession: &Accession) -> Option<&'static str> {
        crate::liftover::aliases::infer_genome_build_from_accession(accession)
    }

    /// Resolve a legacy LOVD gene-model selector (`GENE`, `GENE_v001`) on a
    /// genomic reference to a transcript accession (#500/#637).
    ///
    /// `selector` is the gene-symbol selector verbatim (e.g. `"TIMM8B_v001"`).
    /// `ng_parent` is the versioned `NG_` accession the selector sits on, when
    /// known. When `ng_parent` is `Some` and that exact `NG_` version uniquely
    /// hosts the gene, the hosted transcript takes precedence over the global
    /// reference-standard map (#792); otherwise resolution falls back to that
    /// global map.
    /// Returns the resolved transcript accession (e.g. `"NM_012459.4"`) for a
    /// bare gene or `_v001`, or `None` to decline (an unknown gene, a higher
    /// locus version `_v002…`, or a provider with no gene-model summary
    /// ingested). Declining preserves the input selector rather than emitting a
    /// guessed transcript.
    fn resolve_legacy_gene_selector(
        &self,
        _selector: &str,
        _ng_parent: Option<&Accession>,
    ) -> Option<String> {
        None
    }

    /// Synthesize the transcript selector for a bare-`NG_` `c.` input that
    /// carries no selector, by returning the single transcript the `NG_` parent
    /// hosts when — and only when — that mapping is unambiguous (#923).
    ///
    /// Returns `None` (preserve the bare input) when the parent is absent from
    /// the hosted-transcript map, hosts more than one gene, or its sole gene
    /// hosts more than one transcript. Providers that ingest an
    /// `ng_hosted_transcripts` artifact override this; the default declines.
    fn sole_hosted_transcript(&self, _ng_parent: &Accession) -> Option<String> {
        None
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

/// Blanket implementation for `Arc<T>` where `T: ReferenceProvider + ?Sized`.
///
/// Allows shared-ownership providers (used by `normalize_ferro_parallel`) to
/// be used anywhere a `ReferenceProvider` is expected without unwrapping.
/// This covers both `Arc<ConcreteType>` and `Arc<dyn ReferenceProvider + Send + Sync>`.
impl<T: ReferenceProvider + ?Sized> ReferenceProvider for std::sync::Arc<T> {
    fn get_transcript(&self, id: &str) -> Result<std::sync::Arc<Transcript>, FerroError> {
        (**self).get_transcript(id)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
    ) -> Result<std::sync::Arc<Transcript>, FerroError> {
        (**self).get_transcript_for_variant(variant)
    }

    fn get_transcript_for_accession(
        &self,
        accession: &Accession,
    ) -> Result<std::sync::Arc<Transcript>, FerroError> {
        // Forward so the wrapped provider's parent-context probe survives an
        // Arc wrapper. The default routes through `get_transcript` and would
        // drop the `genomic_context`-aware build resolution.
        (**self).get_transcript_for_accession(accession)
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

    fn genomic_placement(&self, parent: &Accession) -> Option<GenomicPlacement> {
        (**self).genomic_placement(parent)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        // Forward the build hint to the wrapped provider so build-aware lookup
        // survives an Arc/Box wrapper (the default would drop it to the
        // build-agnostic path).
        (**self).genomic_placement_on_build(parent, build)
    }

    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        ng_parent: Option<&Accession>,
    ) -> Option<String> {
        // Forward so legacy gene-model resolution survives an Arc wrapper; the
        // default declines (`None`) and would bypass the wrapped provider.
        (**self).resolve_legacy_gene_selector(selector, ng_parent)
    }

    fn sole_hosted_transcript(&self, ng_parent: &Accession) -> Option<String> {
        // Forward so bare-NG_ selector synthesis survives an Arc wrapper; the
        // default declines (`None`) and would bypass the wrapped provider.
        (**self).sole_hosted_transcript(ng_parent)
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

    fn infer_genome_build(&self, accession: &Accession) -> Option<&'static str> {
        // Forward so the wrapped provider's data-driven table (e.g.
        // MultiFastaProvider's ContigAliases) survives an Arc wrapper; the
        // trait default calls the hardcoded heuristic and would bypass it.
        (**self).infer_genome_build(accession)
    }
}

/// Blanket implementation for boxed trait objects.
///
/// Covers both `Box<dyn ReferenceProvider>` and the wider
/// `Box<dyn ReferenceProvider + Send + Sync>` returned by
/// `create_reference_provider`.
impl<T: ReferenceProvider + ?Sized> ReferenceProvider for Box<T> {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        (**self).get_transcript(id)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        (**self).get_transcript_for_variant(variant)
    }

    fn get_transcript_for_accession(
        &self,
        accession: &Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        (**self).get_transcript_for_accession(accession)
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

    fn genomic_placement(&self, parent: &Accession) -> Option<GenomicPlacement> {
        (**self).genomic_placement(parent)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        // Forward the build hint to the wrapped provider so build-aware lookup
        // survives an Arc/Box wrapper (the default would drop it to the
        // build-agnostic path).
        (**self).genomic_placement_on_build(parent, build)
    }

    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        ng_parent: Option<&Accession>,
    ) -> Option<String> {
        // Forward so legacy gene-model resolution survives a Box wrapper; the
        // default declines (`None`) and would bypass the wrapped provider.
        (**self).resolve_legacy_gene_selector(selector, ng_parent)
    }

    fn sole_hosted_transcript(&self, ng_parent: &Accession) -> Option<String> {
        // Forward so bare-NG_ selector synthesis survives a Box wrapper; the
        // default declines (`None`) and would bypass the wrapped provider.
        (**self).sole_hosted_transcript(ng_parent)
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

    fn infer_genome_build(&self, accession: &Accession) -> Option<&'static str> {
        // Forward so the wrapped provider's data-driven table (e.g.
        // MultiFastaProvider's ContigAliases) survives a Box wrapper; the
        // trait default calls the hardcoded heuristic and would bypass it.
        (**self).infer_genome_build(accession)
    }
}

/// Static assertions: verify that the concrete providers used in production
/// are `Send + Sync` so they can be safely wrapped in `Arc` and shared
/// across threads in `normalize_ferro_parallel`.
///
/// This assertion fires at compile time if any provider adds a non-Send or
/// non-Sync field (e.g. `RefCell`, bare `Rc`, or a raw pointer) in the
/// future.
#[allow(dead_code)]
fn _assert_provider_send_sync() {
    fn assert_send_sync<T: Send + Sync + ?Sized>() {}
    // Trait object form — verifies the trait itself is object-safe and
    // compatible with Send+Sync bounds.
    assert_send_sync::<dyn ReferenceProvider + Send + Sync>();
    // Concrete types used by create_reference_provider in benchmark builds.
    assert_send_sync::<crate::reference::multi_fasta::MultiFastaProvider>();
    assert_send_sync::<crate::reference::mock::MockProvider>();
}

#[cfg(test)]
mod placement_tests {
    use super::*;

    fn placement(strand: Strand) -> GenomicPlacement {
        // Parent base 1 (plus) or the high NC base (minus) anchors at nc 1000;
        // a 9-base span 1000..=1008, so parent runs 1..=9.
        GenomicPlacement {
            nc: Accession::new("NC", "000001", Some(11)),
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1008,
            strand,
        }
    }

    fn placement_on(chrom_version: u32) -> GenomicPlacement {
        GenomicPlacement {
            nc: Accession::new("NC", "000001", Some(chrom_version)),
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1008,
            strand: Strand::Plus,
        }
    }

    #[test]
    fn select_placement_for_build_picks_by_build_and_declines_mismatch() {
        let grch38 = placement_on(11); // NC_000001.11 → GRCh38
        let grch37 = placement_on(10); // NC_000001.10 → GRCh37
        let both = vec![grch38.clone(), grch37.clone()];
        let infer = crate::liftover::aliases::infer_genome_build_from_accession;

        // Explicit build selects that build's placement.
        assert_eq!(
            select_placement_for_build(&both, Some("GRCh38"), infer),
            Some(grch38.clone())
        );
        assert_eq!(
            select_placement_for_build(&both, Some("GRCh37"), infer),
            Some(grch37.clone())
        );
        // No preference → GRCh38 even when listed second.
        assert_eq!(
            select_placement_for_build(&[grch37.clone(), grch38.clone()], None, infer),
            Some(grch38.clone())
        );
        // Build-match guard: an explicit build with no matching placement declines
        // rather than substitute the other build (would mis-anchor).
        let only37 = vec![grch37.clone()];
        assert_eq!(
            select_placement_for_build(&only37, Some("GRCh38"), infer),
            None
        );
        // No preference with only GRCh37 present → return it (never empty if data exists).
        assert_eq!(
            select_placement_for_build(&only37, None, infer),
            Some(grch37)
        );
        assert_eq!(select_placement_for_build(&[], Some("GRCh38"), infer), None);
    }

    #[test]
    fn select_placement_uses_injected_inference() {
        let grch38 = placement_on(11); // NC_000001.11 → GRCh38 (per the hardcoded heuristic)
        let grch37 = placement_on(10); // NC_000001.10 → GRCh37 (per the hardcoded heuristic)
        let both = vec![grch37.clone(), grch38.clone()];
        // Deliberately invert the heuristic: this closure calls NC_000001.10
        // "GRCh38" and everything else `None`. If `select_placement_for_build`
        // ignored the injected closure and reached for the hardcoded helper, it
        // would return `grch38` (the .11 placement) and this assertion would
        // fail — so the test now observes that the injected inferer truly drives
        // selection.
        let infer = |acc: &Accession| {
            if acc.full() == "NC_000001.10" {
                Some("GRCh38")
            } else {
                None
            }
        };
        assert_eq!(
            select_placement_for_build(&both, Some("GRCh38"), infer),
            Some(grch37)
        );
    }

    #[test]
    fn parent_to_nc_plus_strand_is_inverse_of_nc_to_parent() {
        let p = placement(Strand::Plus);
        // Parent base 1 sits at nc 1000; parent base 9 at nc 1008.
        assert_eq!(p.parent_to_nc(1), Some(1000));
        assert_eq!(p.parent_to_nc(4), Some(1003));
        assert_eq!(p.parent_to_nc(9), Some(1008));
        // Round-trips with nc_to_parent across the whole span.
        for nc in p.nc_start..=p.nc_end {
            let parent = p.nc_to_parent(nc).expect("in span");
            assert_eq!(p.parent_to_nc(parent), Some(nc), "round-trip at nc {nc}");
        }
    }

    #[test]
    fn parent_to_nc_minus_strand_reverses_order() {
        let p = placement(Strand::Minus);
        // Minus: parent base 1 anchors at the HIGH nc end (1008).
        assert_eq!(p.parent_to_nc(1), Some(1008));
        assert_eq!(p.parent_to_nc(9), Some(1000));
        for nc in p.nc_start..=p.nc_end {
            let parent = p.nc_to_parent(nc).expect("in span");
            assert_eq!(p.parent_to_nc(parent), Some(nc), "round-trip at nc {nc}");
        }
    }

    #[test]
    fn parent_to_nc_out_of_span_is_none() {
        let p = placement(Strand::Plus);
        // Parent span is 1..=9; 0 and 10 fall outside.
        assert_eq!(p.parent_to_nc(0), None);
        assert_eq!(p.parent_to_nc(10), None);
    }
}

#[cfg(test)]
mod wrapper_forwarding_tests {
    use super::*;

    /// Minimal provider that overrides the methods the `Arc<T>`/`Box<T>`
    /// blanket impls must forward, returning sentinels distinct from each
    /// method's trait default so a test can tell forwarding apart from
    /// fall-through to the default.
    struct OverrideProvider;

    impl ReferenceProvider for OverrideProvider {
        fn get_transcript(&self, _id: &str) -> Result<Arc<Transcript>, FerroError> {
            // The trait default for `get_transcript_for_accession` routes here;
            // a distinct marker lets the test detect that fall-through path.
            Err(FerroError::ReferenceNotFound {
                id: "default-marker".to_string(),
            })
        }

        fn get_sequence(&self, _id: &str, _start: u64, _end: u64) -> Result<String, FerroError> {
            Err(FerroError::ReferenceNotFound {
                id: "unused".to_string(),
            })
        }

        fn get_transcript_for_accession(
            &self,
            _accession: &Accession,
        ) -> Result<Arc<Transcript>, FerroError> {
            Err(FerroError::ReferenceNotFound {
                id: "override-marker".to_string(),
            })
        }

        fn resolve_legacy_gene_selector(
            &self,
            _selector: &str,
            _ng_parent: Option<&Accession>,
        ) -> Option<String> {
            Some("NM_999999.9".to_string())
        }

        fn infer_genome_build(&self, _accession: &Accession) -> Option<&'static str> {
            // Return a sentinel that the hardcoded heuristic would NOT return for
            // an unknown accession (which produces `None`), so a forwarding test
            // can distinguish the override path from trait-default fall-through.
            Some("GRCh38")
        }
    }

    /// Extract the `ReferenceNotFound` id, which carries the sentinel marking
    /// which code path produced the error.
    fn marker(result: Result<Arc<Transcript>, FerroError>) -> String {
        match result {
            Err(FerroError::ReferenceNotFound { id }) => id,
            other => panic!("expected ReferenceNotFound, got {other:?}"),
        }
    }

    #[test]
    fn arc_forwards_overridden_methods() {
        let arc: Arc<dyn ReferenceProvider> = Arc::new(OverrideProvider);
        let accession = Accession::new("NM", "012345", Some(4));
        // Without the forward, the Arc impl falls to the default, which routes
        // through `get_transcript` → "default-marker".
        assert_eq!(
            marker(arc.get_transcript_for_accession(&accession)),
            "override-marker"
        );
        // Without the forward, this drops to the trait default (`None`).
        assert_eq!(
            arc.resolve_legacy_gene_selector("GENE", None),
            Some("NM_999999.9".to_string())
        );
        // Use a bogus NC_ accession the hardcoded table doesn't know → None; the
        // override returns Some("GRCh38"), so forwarding is distinguishable.
        let unknown_nc = Accession::new("NC", "000999", Some(999));
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&unknown_nc),
            None,
            "hardcoded heuristic must return None for the sentinel accession"
        );
        assert_eq!(
            arc.infer_genome_build(&unknown_nc),
            Some("GRCh38"),
            "Arc must forward infer_genome_build to the wrapped provider"
        );
    }

    #[test]
    fn box_forwards_overridden_methods() {
        let boxed: Box<dyn ReferenceProvider> = Box::new(OverrideProvider);
        let accession = Accession::new("NM", "012345", Some(4));
        assert_eq!(
            marker(boxed.get_transcript_for_accession(&accession)),
            "override-marker"
        );
        assert_eq!(
            boxed.resolve_legacy_gene_selector("GENE", None),
            Some("NM_999999.9".to_string())
        );
        // Same forwarding check for Box.
        let unknown_nc = Accession::new("NC", "000999", Some(999));
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&unknown_nc),
            None,
            "hardcoded heuristic must return None for the sentinel accession"
        );
        assert_eq!(
            boxed.infer_genome_build(&unknown_nc),
            Some("GRCh38"),
            "Box must forward infer_genome_build to the wrapped provider"
        );
    }
}
