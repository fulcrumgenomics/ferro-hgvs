//! `VariantProjector` orchestrator.

use crate::data::cdot::CdotTranscript;
use crate::data::mapping::MappingInfo;
use crate::data::projection::Projector;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::{CdsInterval, TxInterval};
use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
use crate::hgvs::variant::{
    is_frameshift, Accession, AlleleVariant, CdsVariant, HgvsVariant, LocEdit, TxVariant,
};
use crate::normalize::{NormalizeConfig, Normalizer};
use crate::project::accession::parse_accession;
use crate::project::edit::transform_edit_for_strand;
use crate::project::protein::{
    affects_initiation_codon, build_initiator_unknown, cds_has_valid_start, predict_indel_protein,
    predict_substitution_protein, read_cds_start_codon, RefProteinBundle,
};
use crate::project::result::VariantProjection;
use crate::reference::transcript::Transcript;
use crate::reference::ReferenceProvider;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

/// Cache key for [`VariantProjector`]'s transcript cache:
/// `(transcript_id, parent_accession)`. The parent component is `None` for
/// non-variant-aware lookups and for inputs without a `genomic_context`; it
/// is `Some` only when the variant-aware path observes an NG/NC parent that
/// could steer the lookup to a non-primary cdot build (issue #332).
type TranscriptCacheKey = (String, Option<String>);

pub struct VariantProjector<P: ReferenceProvider + Clone> {
    projector: Projector,
    provider: P,
    normalizer: Normalizer<P>,
    /// Cache of fetched transcripts keyed by `(transcript_id, parent_accession)`.
    ///
    /// `project_single_inner` looks up each overlapping transcript via the
    /// provider once per variant; for `project_all` workloads (hundreds of
    /// overlapping transcripts per variant on dense chromosomes) the same
    /// transcript ID is fetched repeatedly. Caching collapses N→1 fetches
    /// per unique key over the projector's lifetime.
    ///
    /// The parent-accession component is non-`None` only on the variant-aware
    /// lookup path (`cached_get_transcript_for_variant`), where the variant's
    /// `genomic_context` (NG/NC parent) can steer
    /// `ReferenceProvider::get_transcript_for_variant` to a different cdot
    /// build for the same `transcript_id` (issue #332). Inputs without a
    /// `genomic_context` and the legacy `cached_get_transcript` entry point
    /// both cache under `(transcript_id, None)`, so they share entries.
    transcript_cache: RwLock<HashMap<TranscriptCacheKey, Arc<Transcript>>>,
    /// Cache of the translated reference protein (ref CDS + ref protein +
    /// ref-with-stop) keyed by the same `(transcript_id, parent_accession)`
    /// identity used by `transcript_cache`. Each entry is stable for the life
    /// of the projector and is consumed by `predict_indel_protein` to avoid
    /// retranslating the full CDS on every indel.
    ///
    /// The parent-aware key matters because one `transcript_id` can resolve to
    /// different transcript records (and therefore different reference CDS
    /// bases) under different NG/NC parents (issue #332). Keying by
    /// `transcript_id` alone would let the first resolved build poison indel
    /// protein predictions for projections on the other build.
    ref_protein_cache: RwLock<HashMap<TranscriptCacheKey, Arc<RefProteinBundle>>>,
}

impl<P: ReferenceProvider + Clone> VariantProjector<P> {
    pub fn new(projector: Projector, provider: P) -> Self {
        let normalizer = Normalizer::with_config(provider.clone(), NormalizeConfig::default());
        Self {
            projector,
            provider,
            normalizer,
            transcript_cache: RwLock::new(HashMap::new()),
            ref_protein_cache: RwLock::new(HashMap::new()),
        }
    }

    /// Non-variant-aware transcript fetch with caching. Kept for inline tests
    /// that need to populate a `Transcript` for `cached_ref_translation`
    /// without constructing an `HgvsVariant`. Hot-path callers go through
    /// [`Self::cached_get_transcript_for_variant`] so an NG/NC-parented input
    /// resolves to the correct cdot build (issue #332).
    #[cfg(test)]
    fn cached_get_transcript(&self, transcript_id: &str) -> Result<Arc<Transcript>, FerroError> {
        let key = (transcript_id.to_string(), None);
        if let Some(tx) = self
            .transcript_cache
            .read()
            .expect("transcript cache poisoned")
            .get(&key)
        {
            return Ok(Arc::clone(tx));
        }
        let tx = self.provider.get_transcript(transcript_id)?;
        let mut guard = self
            .transcript_cache
            .write()
            .expect("transcript cache poisoned");
        // `tx` is already an `Arc<Transcript>` from the provider; insert directly.
        let entry = guard.entry(key).or_insert_with(|| tx);
        Ok(Arc::clone(entry))
    }

    /// Variant-aware transcript lookup with caching.
    ///
    /// Routes through [`ReferenceProvider::get_transcript_for_variant`] so an
    /// NG/NC-parented input picks the build-correct cdot transcript (issue
    /// #332). The cache key is `(transcript_id, parent_accession)` because
    /// the same `transcript_id` can resolve to different transcripts under
    /// different parents.
    ///
    /// On `ReferenceNotFound` from the variant-aware call, falls back to the
    /// bare provider lookup (matches the prior `tx_for_codon_with_fallback`
    /// contract). Other provider errors propagate.
    /// Build the parent-aware cache key component for a variant.
    ///
    /// Returns `Some(parent.to_string())` when the variant's accession carries
    /// a `genomic_context` (NG/NC parent), `None` otherwise. Shared by
    /// [`Self::cached_get_transcript_for_variant`] and
    /// [`Self::cached_ref_translation`] so both caches use byte-identical keys
    /// for the same variant.
    fn parent_key_for(variant: &HgvsVariant) -> Option<String> {
        variant
            .accession()
            .and_then(|a| a.genomic_context.as_deref())
            .map(|gc| gc.to_string())
    }

    /// Reject parentless c./n./r. fan-out inputs up front.
    ///
    /// `project_normalized_all` seeds a stab query against the contig
    /// derived from cdot's exon table, then projects through every
    /// overlapping transcript via `project_single_inner`. For c./n./r.
    /// inputs `project_single_inner` calls `project_to_genomic`, which
    /// requires an explicit NG/NC parent on `accession.genomic_context`
    /// (see #327). Without the parent every per-transcript projection
    /// fails with `UnsupportedProjection { "no parent reference..." }`,
    /// and the fan-out loop's `log::trace!`-and-drop converts the user
    /// error into `Ok([])` — silently turning a missing-context error
    /// into "no overlaps." Surface the error here, before we touch the
    /// stab-query path (#389 follow-up).
    fn require_parent_for_fanout(
        accession: &crate::hgvs::variant::Accession,
        axis: &str,
    ) -> Result<(), FerroError> {
        if accession.genomic_context.is_none() {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "project_all requires a parent reference (genomic_context) on {} \
                     inputs to resolve back to g.; accession {} has none",
                    axis,
                    accession.full()
                ),
            });
        }
        Ok(())
    }

    /// Infer the genome build to consult cdot under for `variant`.
    ///
    /// Examines the variant's own accession first: a g. variant on
    /// `NC_*.10` carries GRCh37, on `NC_*.11` carries GRCh38, and an
    /// assembly-style `GRCh37(chr1):g.…` carries GRCh37 directly. If the
    /// variant has no inherent build (e.g. a c. variant on a bare NM)
    /// but does have a `genomic_context` (NG/NC parent), the parent is
    /// consulted. Anything else returns `None` so the caller falls back
    /// to the build-agnostic primary-cdot path (issue #389).
    fn build_hint_for_variant(variant: &HgvsVariant) -> Option<&'static str> {
        use crate::liftover::aliases::infer_genome_build_from_accession;
        let acc = variant.accession()?;
        if let Some(b) = infer_genome_build_from_accession(acc) {
            return Some(b);
        }
        if let Some(parent) = acc.genomic_context.as_deref() {
            return infer_genome_build_from_accession(parent);
        }
        None
    }

    /// Build-aware cdot transcript lookup shared by every projection
    /// entry point: routes through
    /// [`CdotMapper::get_transcript_on_build`] when `build_hint` is
    /// `Some`, falling back to the primary-build
    /// [`CdotMapper::get_transcript`] otherwise. Returns
    /// [`FerroError::ReferenceNotFound`] if the requested build's view
    /// is absent from cdot (#389).
    fn cdot_tx_with_build_hint(
        &self,
        transcript_id: &str,
        build_hint: Option<&'static str>,
    ) -> Result<&CdotTranscript, FerroError> {
        match build_hint {
            Some(b) => self
                .projector
                .mapper()
                .cdot()
                .get_transcript_on_build(transcript_id, b),
            None => self.projector.mapper().cdot().get_transcript(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })
    }

    /// Extract the contig name and a representative 1-based genomic
    /// position from an already-normalized [`HgvsVariant`], for use as
    /// the input to `Projector::project`'s stab query.
    ///
    /// For [`HgvsVariant::Allele`] the first inner variant is used.
    ///
    /// # Per-axis behavior
    ///
    /// - **Genome**: contig is taken from the variant's accession
    ///   (`chromosome` for assembly-style refs, `full()` otherwise);
    ///   position is the start of the genomic interval.
    /// - **Cds / Tx / Rna**: cdot is consulted to resolve the
    ///   transcript's contig and to project the start position to
    ///   genome — `cds_to_genome_on_build` for c. (handles
    ///   intronic offsets, 5'UTR negative bases, and 3'UTR `utr3`
    ///   natively); `tx_to_genome` for n./r. simple exonic positions.
    ///   For n./r. positions that `tx_to_genome` cannot resolve
    ///   precisely (intronic offsets, downstream/utr3 markers,
    ///   non-positive bases) the function falls back to the
    ///   transcript's first-exon genomic start — the stab query still
    ///   lands inside the transcript, and the per-transcript
    ///   projection downstream re-derives the accurate coordinate. The
    ///   contig is taken from `cdot_tx.contig` (production cdot
    ///   typically keys this by the RefSeq genomic accession, with
    ///   `chr*` aliases handled by `resolve_contig` inside the stab
    ///   query).
    ///
    /// Removes the g.-only gate that pre-#389 defeated PR #379's
    /// fan-out widening. Build hints are inferred from the input's
    /// `genomic_context` parent (or assembly tag) the same way as the
    /// rest of the projector.
    fn extract_contig_and_pos(&self, variant: &HgvsVariant) -> Result<(String, u64), FerroError> {
        let effective = match variant {
            HgvsVariant::Allele(allele) => {
                allele
                    .variants
                    .first()
                    .ok_or_else(|| FerroError::UnsupportedProjection {
                        reason: "cannot project an empty allele to all transcripts".to_string(),
                    })?
            }
            other => other,
        };

        match effective {
            HgvsVariant::Genome(gv) => {
                let contig = if let Some(chr) = &gv.accession.chromosome {
                    chr.to_string()
                } else {
                    gv.accession.full()
                };
                let pos = gv
                    .loc_edit
                    .location
                    .start
                    .inner()
                    .cloned()
                    .ok_or_else(|| FerroError::InvalidCoordinates {
                        msg: "genomic interval start is unknown".to_string(),
                    })?
                    .base;
                Ok((contig, pos))
            }
            HgvsVariant::Cds(c) => {
                Self::require_parent_for_fanout(&c.accession, "c.")?;
                let transcript_id = c.accession.transcript_accession();
                let build_hint = Self::build_hint_for_variant(effective);
                let cdot_tx = self.cdot_tx_with_build_hint(&transcript_id, build_hint)?;
                let start_cds = c.loc_edit.location.start.inner().cloned().ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: "CDS interval start is unknown".to_string(),
                    }
                })?;
                let pos = self
                    .projector
                    .mapper()
                    .cds_to_genome_on_build(&transcript_id, &start_cds, build_hint)?
                    .variant
                    .base;
                Ok((cdot_tx.contig.clone(), pos))
            }
            HgvsVariant::Tx(t) => {
                Self::require_parent_for_fanout(&t.accession, "n.")?;
                let transcript_id = t.accession.transcript_accession();
                let build_hint = Self::build_hint_for_variant(effective);
                let cdot_tx = self.cdot_tx_with_build_hint(&transcript_id, build_hint)?;
                let start_tx = *t.loc_edit.location.start.inner().ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: "Tx interval start is unknown".to_string(),
                    }
                })?;
                let pos =
                    nr_representative_genome_pos(cdot_tx, start_tx.base, start_tx.offset, false)?;
                Ok((cdot_tx.contig.clone(), pos))
            }
            HgvsVariant::Rna(r) => {
                Self::require_parent_for_fanout(&r.accession, "r.")?;
                let transcript_id = r.accession.transcript_accession();
                let build_hint = Self::build_hint_for_variant(effective);
                let cdot_tx = self.cdot_tx_with_build_hint(&transcript_id, build_hint)?;
                let start_rna = *r.loc_edit.location.start.inner().ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: "RNA interval start is unknown".to_string(),
                    }
                })?;
                let pos = nr_representative_genome_pos(
                    cdot_tx,
                    start_rna.base,
                    start_rna.offset,
                    start_rna.utr3,
                )?;
                Ok((cdot_tx.contig.clone(), pos))
            }
            HgvsVariant::Protein(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept protein (p.) inputs".to_string(),
            }),
            HgvsVariant::Mt(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept mitochondrial (m.) inputs".to_string(),
            }),
            HgvsVariant::Circular(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept circular (o.) inputs".to_string(),
            }),
            HgvsVariant::RnaFusion(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept RNA fusion inputs".to_string(),
            }),
            HgvsVariant::GenomeRing(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept genome ring inputs".to_string(),
            }),
            HgvsVariant::Supernumerary(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept supernumerary (sup) inputs".to_string(),
            }),
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => {
                Err(FerroError::UnsupportedProjection {
                    reason: "project_all does not accept null/unknown allele inputs".to_string(),
                })
            }
            // Allele was already unwrapped to its first inner variant above;
            // if we reach this arm with one it's because the inner variant
            // didn't match any of the supported axes (e.g. p./m./o.).
            HgvsVariant::Allele(_) => Err(FerroError::UnsupportedProjection {
                reason:
                    "project_all does not accept alleles whose first inner variant is not g./c./n./r."
                        .to_string(),
            }),
        }
    }

    fn cached_get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<Arc<Transcript>, FerroError> {
        let parent_key = Self::parent_key_for(variant);
        let key = (transcript_id.to_string(), parent_key);
        if let Some(tx) = self
            .transcript_cache
            .read()
            .expect("transcript cache poisoned")
            .get(&key)
        {
            return Ok(Arc::clone(tx));
        }
        let tx = match self.provider.get_transcript_for_variant(variant) {
            Ok(tx) => tx,
            Err(FerroError::ReferenceNotFound { .. }) => {
                self.provider.get_transcript(transcript_id)?
            }
            Err(e) => return Err(e),
        };
        let mut guard = self
            .transcript_cache
            .write()
            .expect("transcript cache poisoned");
        // `tx` is already an `Arc<Transcript>` from the provider; insert directly.
        let entry = guard.entry(key).or_insert_with(|| tx);
        Ok(Arc::clone(entry))
    }

    /// Build (or retrieve) the cached `RefProteinBundle` for a transcript.
    ///
    /// The reference CDS and its translation are stable for a given
    /// `(transcript_id, parent_accession)` pair — `predict_indel_protein`
    /// would otherwise recompute them on every variant. Caching collapses
    /// that to one translation per (id, parent) over the projector's
    /// lifetime.
    ///
    /// The cache is keyed by the same parent-aware identity as
    /// [`Self::cached_get_transcript_for_variant`] so that an NG/NC-parented
    /// input does not reuse the bundle built for a different cdot build of
    /// the same `transcript_id` (issue #332).
    pub(crate) fn cached_ref_translation(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
        transcript: &Transcript,
    ) -> Result<Arc<RefProteinBundle>, FerroError> {
        let parent_key = Self::parent_key_for(variant);
        let key = (transcript_id.to_string(), parent_key);
        if let Some(b) = self
            .ref_protein_cache
            .read()
            .expect("ref-protein cache poisoned")
            .get(&key)
        {
            return Ok(Arc::clone(b));
        }
        let bundle = RefProteinBundle::from_transcript(transcript)?;
        let mut guard = self
            .ref_protein_cache
            .write()
            .expect("ref-protein cache poisoned");
        let entry = guard.entry(key).or_insert_with(|| Arc::new(bundle));
        Ok(Arc::clone(entry))
    }

    pub fn with_normalize_config(mut self, config: NormalizeConfig) -> Self {
        self.normalizer = Normalizer::with_config(self.provider.clone(), config);
        self
    }

    /// Parse, normalize, and project an HGVS string onto a transcript.
    pub fn project(
        &self,
        hgvs_string: &str,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        let variant = crate::parse_hgvs(hgvs_string)?;
        self.project_variant(&variant, transcript_id)
    }

    /// Normalize and project an already-parsed g. variant onto a transcript.
    ///
    /// The variant is normalized first; for pre-normalized variants use
    /// [`project_normalized`] to skip the redundant normalization step.
    pub fn project_variant(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // 1. Normalize the genomic variant. The normalizer is built once at
        // construction time so we don't clone the (potentially heavy) provider
        // on every call.
        let normalized = self.normalizer.normalize(variant)?;
        self.project_variant_inner(&normalized, transcript_id)
    }

    /// Project an already-normalized g. variant onto a transcript, skipping the
    /// normalization step.
    ///
    /// Callers that pre-normalize once and then project against many transcripts
    /// should use this method to avoid the cost of re-normalization.
    ///
    /// **Warning**: passing a non-normalized variant will produce coordinates
    /// that are technically valid but may not match other tools' canonical form.
    pub fn project_normalized(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        self.project_variant_inner(variant, transcript_id)
    }

    /// Parse, normalize, and project an HGVS string onto ALL overlapping
    /// transcripts, returning results in clinical priority order (MANE Select
    /// first, then Plus Clinical, then canonical, then longest CDS).
    ///
    /// Returns an empty `Vec` when the variant overlaps no known transcripts.
    /// Individual transcript errors are logged at trace level and silently
    /// skipped so that a single bad transcript does not abort the whole call.
    pub fn project_all(&self, hgvs_string: &str) -> Result<Vec<VariantProjection>, FerroError> {
        let variant = crate::parse_hgvs(hgvs_string)?;
        self.project_variant_all(&variant)
    }

    /// Normalize and project an already-parsed g. variant onto ALL overlapping
    /// transcripts.
    ///
    /// See [`project_all`] for ordering and error-handling semantics.
    pub fn project_variant_all(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VariantProjection>, FerroError> {
        // 1. Normalize once via the cached normalizer (built at construction time).
        let normalized = self.normalizer.normalize(variant)?;
        self.project_normalized_all(&normalized)
    }

    /// Derive the genomic LRG parent (`LRG_<n>`) of an LRG transcript or protein
    /// accession (`LRG_<n>t<m>` / `LRG_<n>p<m>`).
    ///
    /// LRG accessions embed the genomic number plus a `t`/`p` suffix, so the
    /// genomic reference is the leading digits — a structural fact of the
    /// accession, not a cdot inference. Returns `None` for non-LRG accessions
    /// and for a bare genomic LRG (all digits, no suffix), which is already its
    /// own genome reference.
    fn lrg_genomic_parent(accession: &Accession) -> Option<Accession> {
        if !accession.is_lrg() {
            return None;
        }
        let number = &*accession.number;
        let genomic_digits: String = number.chars().take_while(|c| c.is_ascii_digit()).collect();
        if genomic_digits.is_empty() || genomic_digits.len() == number.len() {
            return None;
        }
        Some(Accession::new("LRG", genomic_digits, None))
    }

    /// Project a transcript-coordinate variant (`c.`/`n.`/`r.`) onto its parent
    /// genomic reference and return an [`HgvsVariant::Genome`].
    ///
    /// The output uses the parent NG/NC accession stored in the input's
    /// `Accession.genomic_context`. The function is idempotent on `Genome`
    /// input (returned unchanged).
    ///
    /// # Errors
    ///
    /// Returns [`FerroError::UnsupportedProjection`] for:
    /// - `p.` / `m.` / `o.` / RNA-fusion / `Allele` / `NullAllele` / `UnknownAllele`
    ///   inputs (see #328 for allele support),
    /// - transcript-coordinate inputs whose `Accession.genomic_context` is
    ///   absent (no parent NG/NC reference),
    /// - `?` position sentinels in the start or end coordinate.
    ///
    /// Returns [`FerroError::InvalidCoordinates`] when an endpoint is a
    /// compound `(a_b)` `UncertainBoundary::Range` (no single position).
    ///
    /// Returns [`FerroError::ReferenceNotFound`] when the named transcript
    /// is absent from the cdot mapper.
    ///
    /// The output ref-base is **not** validated against the genomic reference;
    /// mismatches surface downstream in `Normalizer`.
    ///
    /// **Limitation (`r.`):** plus-strand RNA inputs whose edit carries
    /// `Base::U` are forwarded unchanged into the g. output, producing an
    /// invalid DNA emission. Callers should pre-translate U→T before
    /// constructing the r. variant. Tracked separately.
    pub fn project_to_genomic(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
        use crate::hgvs::variant::{GenomeVariant, HgvsVariant, LocEdit};
        use crate::reference::Strand as RefStrand;

        // The `UncertainBoundary<T>` → `T` resolver lives at module scope
        // (`resolve_uncertain_boundary`) so the direct bare-NM_ path shares the
        // exact `?`-vs-range classification used here.

        // 0. Reject variant kinds that have no genomic coordinate counterpart
        //    in this iteration (see #328 for allele support).
        match variant {
            HgvsVariant::Genome(_) => return Ok(variant.clone()),
            HgvsVariant::Protein(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support protein (p.) inputs".to_string(),
                })
            }
            HgvsVariant::Mt(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support mitochondrial (m.) inputs"
                        .to_string(),
                })
            }
            HgvsVariant::Circular(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support circular (o.) inputs".to_string(),
                })
            }
            HgvsVariant::RnaFusion(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support RNA fusion inputs".to_string(),
                })
            }
            HgvsVariant::GenomeRing(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support genome ring inputs".to_string(),
                })
            }
            HgvsVariant::Supernumerary(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support supernumerary (sup) inputs"
                        .to_string(),
                })
            }
            HgvsVariant::Allele(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not yet support allele inputs; see #328"
                        .to_string(),
                })
            }
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support null/unknown allele markers"
                        .to_string(),
                })
            }
            HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_) => {}
        }

        // 1. Pull the common pieces (accession, gene_symbol, edit, raw start/end
        //    as CdsPos) out of each transcript-coord variant kind. RNA and Tx
        //    positions have the same shape as CdsPos (base / offset / utr3-or-
        //    downstream), so we normalize to CdsPos here. Coding positions
        //    route through `data::mapping::CoordinateMapper::cds_to_genome`
        //    (which also handles intronic offsets and the `utr3` flag); the
        //    non-coding path uses `CdotTranscript::tx_to_genome` directly.
        let (accession, gene_symbol, edit_mu, start_cds, end_cds) = match variant {
            HgvsVariant::Cds(v) => {
                let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "c.", "start")?;
                let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "c.", "end")?;
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    s,
                    e,
                )
            }
            HgvsVariant::Tx(v) => {
                let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "n.", "start")?;
                let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "n.", "end")?;
                let to_cds = |p: TxPos| CdsPos {
                    base: p.base,
                    offset: p.offset,
                    utr3: p.downstream,
                    special: None,
                };
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    to_cds(s),
                    to_cds(e),
                )
            }
            HgvsVariant::Rna(v) => {
                let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "r.", "start")?;
                let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "r.", "end")?;
                // RnaPos shape mirrors CdsPos exactly.
                let to_cds = |p: crate::hgvs::location::RnaPos| CdsPos {
                    base: p.base,
                    offset: p.offset,
                    utr3: p.utr3,
                    special: None,
                };
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    to_cds(s),
                    to_cds(e),
                )
            }
            _ => unreachable!("variant kind already filtered above"),
        };

        // 2. Reject `?` position sentinels explicitly — these have base == 0
        //    and no offset, which would otherwise propagate into the genomic
        //    mapper as a meaningless coordinate.
        if start_cds.is_unknown()
            || start_cds.is_special()
            || end_cds.is_unknown()
            || end_cds.is_special()
        {
            return Err(FerroError::UnsupportedProjection {
                reason: "project_to_genomic cannot resolve `?` or special (pter/qter/cen) \
                         position sentinels"
                    .to_string(),
            });
        }

        // 3. Resolve the genomic parent reference. Normally this is an explicit
        //    NG/NC parent on `genomic_context` — per #327 we do NOT synthesize a
        //    parent from cdot. The one exception is an LRG transcript/protein
        //    input (`LRG_<n>t<m>` / `LRG_<n>p<m>`): its genomic parent `LRG_<n>`
        //    is determined *structurally* by the accession itself, not inferred
        //    from cdot, so we derive it when no explicit parent is given (#480).
        let parent = match accession.genomic_context.as_deref().cloned() {
            Some(p) => p,
            None => Self::lrg_genomic_parent(&accession).ok_or_else(|| {
                FerroError::UnsupportedProjection {
                    reason: format!(
                        "input variant has no parent reference (genomic_context) on accession {}; \
                         project_to_genomic requires an explicit NG/NC parent (see #327)",
                        accession.full()
                    ),
                }
            })?,
        };

        // 4. Resolve the transcript via the cdot mapper (the same backing
        //    store used by the g. → c./n. path).  Working directly off cdot
        //    keeps us symmetric with `project_single_inner` and avoids a
        //    second `provider.get_transcript` load whose exon entries may
        //    lack genomic coordinates. When the NG/NC parent's version
        //    disambiguates a genome build (NC_*.10 → GRCh37, NC_*.11 →
        //    GRCh38) we ask cdot for that build's view so multi-build cdot
        //    loads project against the correct alignment (issue #389).
        let transcript_id = accession.transcript_accession();
        let build_hint = crate::liftover::aliases::infer_genome_build_from_accession(&parent);
        let cdot_mapper = self.projector.mapper();
        let cdot_tx = match build_hint {
            Some(b) => cdot_mapper
                .cdot()
                .get_transcript_on_build(&transcript_id, b),
            None => cdot_mapper.cdot().get_transcript(&transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;
        let strand = cdot_tx.strand;
        if strand == RefStrand::Unknown {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "transcript {} has unknown strand; cannot project to g.",
                    transcript_id
                ),
            });
        }

        // 5. Pick the right primitive: for coding transcripts use the
        //    CDS-aware path from `data::mapping` (handles intronic offsets,
        //    UTR markers, and strand swap). For non-coding transcripts the
        //    input `base` is the 1-based tx_pos directly; convert it to the
        //    0-based form cdot exposes and use `CdotTranscript::tx_to_genome`.
        let is_coding = cdot_tx.cds_start.is_some() && cdot_tx.cds_end.is_some();
        let map_pos = |p: CdsPos| -> Result<u64, FerroError> {
            if is_coding {
                let result = cdot_mapper.cds_to_genome_on_build(&transcript_id, &p, build_hint)?;
                Ok(result.variant.base)
            } else {
                if p.offset.is_some() {
                    // Intronic n. offsets require coding-style exon-boundary
                    // arithmetic; we don't currently expose that primitive
                    // for non-coding transcripts. Surface a clear error.
                    return Err(FerroError::UnsupportedProjection {
                        reason: format!(
                            "project_to_genomic does not yet support intronic offsets on \
                             non-coding transcripts ({}); see #332",
                            transcript_id
                        ),
                    });
                }
                // `base` is 1-based on TxPos / RnaPos; cdot's tx coords are
                // 0-based, so subtract 1.
                if p.base <= 0 {
                    return Err(FerroError::InvalidCoordinates {
                        msg: format!(
                            "non-coding transcript {} position must be positive, got {}",
                            transcript_id, p.base
                        ),
                    });
                }
                let tx_pos_0based = (p.base - 1) as u64;
                cdot_tx
                    .tx_to_genome(tx_pos_0based)
                    .ok_or_else(|| FerroError::InvalidCoordinates {
                        msg: format!(
                            "tx position {} not in any exon of {}",
                            p.base, transcript_id
                        ),
                    })
            }
        };

        let start_g = map_pos(start_cds)?;
        let end_g = map_pos(end_cds)?;

        // 6. Build the genomic interval. On minus strand the c./n./r. interval
        //    runs anti-parallel to the genome, so the smaller genomic coord
        //    becomes the start.
        let (g_start_pos, g_end_pos) = match strand {
            RefStrand::Plus => (GenomePos::new(start_g), GenomePos::new(end_g)),
            RefStrand::Minus => (GenomePos::new(end_g), GenomePos::new(start_g)),
            RefStrand::Unknown => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "transcript {} has unknown strand; cannot project to g.",
                        transcript_id
                    ),
                });
            }
        };
        let g_interval = GenomeInterval::new(g_start_pos, g_end_pos);

        // 7. Transform the edit for the transcript strand. The c./n./r. edit
        //    reads on the transcript's sense strand; on minus strand we need
        //    to reverse-complement to put it back onto the genomic reference.
        //    We extract the concrete `NaEdit` (preserving `Mu` certainty) and
        //    apply `transform_edit_for_strand` from the project::edit module.
        let edit_inner =
            edit_mu
                .inner()
                .cloned()
                .ok_or_else(|| FerroError::UnsupportedProjection {
                    reason: "project_to_genomic requires a concrete edit (not `?`)".to_string(),
                })?;
        let g_edit_inner: NaEdit = transform_edit_for_strand(&edit_inner, strand);
        // Preserve the Mu certainty (Certain vs Uncertain).
        let g_edit_mu = edit_mu.map(|_| g_edit_inner);

        // 8. Re-anchor into the parent's own frame when the parent is an
        //    NG_/LRG_ whose chromosomal placement is known (#480). cdot
        //    resolved the coordinates on the chromosome (`NC_`), but `parent`
        //    is an NG_/LRG_ whose own frame differs: apply the affine
        //    NC→parent transform, and reverse-complement the edit when the
        //    parent runs antiparallel to the chromosome. Without a placement
        //    (e.g. an `NC_` parent, or an `NG_` whose placement has not been
        //    ingested) the chromosome coordinates are kept as-is — the prior
        //    behavior — rather than failing the projection.
        let (g_interval, g_edit_mu) = match self.provider.genomic_placement(&parent) {
            Some(placement) => {
                // `g_start_pos.base <= g_end_pos.base` by construction above.
                let nc_lo = g_start_pos.base;
                let nc_hi = g_end_pos.base;
                match (placement.nc_to_parent(nc_lo), placement.nc_to_parent(nc_hi)) {
                    (Some(p_a), Some(p_b)) => {
                        // A minus-strand placement reverses coordinate order.
                        let (p_lo, p_hi) = if p_a <= p_b { (p_a, p_b) } else { (p_b, p_a) };
                        let parent_interval =
                            GenomeInterval::new(GenomePos::new(p_lo), GenomePos::new(p_hi));
                        let parent_edit_mu = if placement.strand == RefStrand::Minus {
                            g_edit_mu
                                .map(|inner| transform_edit_for_strand(&inner, RefStrand::Minus))
                        } else {
                            g_edit_mu
                        };
                        (parent_interval, parent_edit_mu)
                    }
                    // An endpoint fell outside the placed span — keep the
                    // chromosome coordinates rather than emit an out-of-frame
                    // position.
                    _ => (g_interval, g_edit_mu),
                }
            }
            None => (g_interval, g_edit_mu),
        };

        let g_variant = GenomeVariant {
            accession: parent,
            gene_symbol,
            loc_edit: LocEdit::with_uncertainty(g_interval, g_edit_mu),
        };
        Ok(HgvsVariant::Genome(g_variant))
    }

    /// Project an already-normalized g. variant onto ALL overlapping
    /// transcripts, skipping re-normalization.
    ///
    /// Callers that pre-normalize once and then fan-out across transcripts
    /// should use this method.
    pub fn project_normalized_all(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VariantProjection>, FerroError> {
        // 2. Extract contig + first-base position from the normalized variant.
        //    For Allele variants we use the first inner variant's accession.
        //    c./n./r. inputs are resolved via cdot (contig from the
        //    transcript's exon table, position from the start-of-range
        //    transcript→genome projection) so PR #379's per-transcript
        //    widening reaches the fan-out path here (issue #389 item 3).
        let (contig, pos) = self.extract_contig_and_pos(variant)?;

        // 3. Find overlapping transcripts via the Projector (sorted by priority).
        let projection_result = self.projector.project(&contig, pos)?;

        // 4. Project against each overlapping transcript.
        let mut results = Vec::with_capacity(projection_result.projections.len());
        for tx_proj in &projection_result.projections {
            match self.project_variant_inner(variant, &tx_proj.transcript_id) {
                Ok(vp) => results.push(vp),
                Err(e) => {
                    // Skip transcripts that fail — log at trace level only.
                    log::trace!(
                        "project_normalized_all: skipping {} for {}: {}",
                        tx_proj.transcript_id,
                        variant,
                        e
                    );
                }
            }
        }

        Ok(results)
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /// Core projection logic, operating on a variant that is assumed to be
    /// already normalized.  Does NOT re-normalize.
    fn project_variant_inner(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Dispatch on variant kind.
        match variant {
            HgvsVariant::Allele(allele) => {
                self.project_allele_inner(allele, variant, transcript_id)
            }
            _ => self.project_single_inner(variant, transcript_id),
        }
    }

    /// Project a compound [`AlleleVariant`] by recursively projecting each
    /// inner variant and combining the results.
    fn project_allele_inner(
        &self,
        allele: &AlleleVariant,
        original: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Empty allele: pass through with no coding/protein, but validate the
        // transcript ID — every other path through `project_variant_inner` errors
        // out for unknown transcripts (via `project_single_inner`), and this
        // branch must behave the same to avoid silently returning bogus
        // projections for typo'd or missing accessions.
        if allele.variants.is_empty() {
            // Honor the input's parent NG/NC (or assembly tag) when
            // resolving the gene_name from cdot, so an NG/NC-parented
            // allele on a multi-build cdot picks the correct build's
            // transcript record (issue #389).
            let build_hint = Self::build_hint_for_variant(original);
            let gene_symbol = match build_hint {
                Some(b) => self
                    .projector
                    .mapper()
                    .cdot()
                    .get_transcript_on_build(transcript_id, b),
                None => self.projector.mapper().cdot().get_transcript(transcript_id),
            }
            .ok_or_else(|| FerroError::ReferenceNotFound {
                id: transcript_id.to_string(),
            })?
            .gene_name
            .clone();
            return Ok(VariantProjection {
                // No inner variant to derive a genomic form from; report `None`
                // rather than wrapping the (possibly non-genomic) input allele
                // into `.genomic`. Consistent with `coding`/`protein`, which are
                // already `None` for the empty allele (#508 review).
                genomic: None,
                coding: None,
                // Empty allele: no inner member to derive an n. form from.
                noncoding: None,
                protein: None,
                transcript_id: transcript_id.to_string(),
                gene_symbol,
                is_frameshift: false,
                is_intronic: false,
                is_utr: false,
            });
        }

        // Recursively project each inner variant.
        let mut inner_projections = Vec::with_capacity(allele.variants.len());
        for inner in &allele.variants {
            let proj = self.project_variant_inner(inner, transcript_id)?;
            inner_projections.push(proj);
        }

        // Aggregate flags.
        let is_frameshift = inner_projections.iter().any(|p| p.is_frameshift);
        let is_intronic = inner_projections.iter().any(|p| p.is_intronic);
        let is_utr = inner_projections.iter().any(|p| p.is_utr);

        // gene_symbol from any projection that has one.
        let gene_symbol = inner_projections.iter().find_map(|p| p.gene_symbol.clone());

        // Build the coding allele from inner c./n. variants.
        let coding_variants: Option<Vec<HgvsVariant>> =
            inner_projections.iter().map(|p| p.coding.clone()).collect();
        let coding = coding_variants
            .map(|variants| HgvsVariant::Allele(AlleleVariant::new(variants, allele.phase)));

        // Build the n. (transcript-relative) allele the same all-or-nothing way:
        // only when EVERY member carries an n. form (mirrors the coding rule).
        let noncoding_variants: Option<Vec<HgvsVariant>> = inner_projections
            .iter()
            .map(|p| p.noncoding.clone())
            .collect();
        let noncoding = noncoding_variants
            .map(|variants| HgvsVariant::Allele(AlleleVariant::new(variants, allele.phase)));

        // Build the protein allele only if ALL inner projections have a protein.
        let all_have_protein = inner_projections.iter().all(|p| p.protein.is_some());
        let protein = if all_have_protein {
            let protein_variants: Vec<HgvsVariant> = inner_projections
                .iter()
                .filter_map(|p| p.protein.clone())
                .collect();
            Some(HgvsVariant::Allele(AlleleVariant::new(
                protein_variants,
                allele.phase,
            )))
        } else {
            None
        };

        // Derive `.genomic` from the inner projections, not the raw `original`:
        // for a c./n./r. allele `original` is a non-genomic allele, and a
        // bare-NM_ inner projection has no genomic form at all. Mirror the
        // protein rule — build a g. allele only when EVERY member carries a
        // genomic form; otherwise report `None` so a genome-less projection
        // honestly says so rather than stuffing a transcript-coordinate allele
        // into `.genomic` (#508 review).
        let all_have_genomic = inner_projections.iter().all(|p| p.genomic.is_some());
        let genomic = if all_have_genomic {
            let genomic_variants: Vec<HgvsVariant> = inner_projections
                .iter()
                .filter_map(|p| p.genomic.clone())
                .collect();
            Some(HgvsVariant::Allele(AlleleVariant::new(
                genomic_variants,
                allele.phase,
            )))
        } else {
            None
        };

        Ok(VariantProjection {
            genomic,
            coding,
            noncoding,
            protein,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift,
            is_intronic,
            is_utr,
        })
    }

    /// Project a single (non-allele) variant, assuming it has already been
    /// normalized.
    ///
    /// Accepts g./c./n./r. inputs. Non-genomic inputs are first
    /// projected back to their parent NG/NC reference via
    /// [`Self::project_to_genomic`], then re-enter the g.-anchored
    /// projection path. The parent reference must be present in the
    /// input's `Accession.genomic_context` (per #327); inputs without
    /// it surface a clear "no parent reference" diagnostic.
    ///
    /// **Transcript-id constraint.** For c./n./r. inputs the caller's
    /// `transcript_id` must equal the input's own
    /// `accession.transcript_accession()`. Mixing the two (e.g. passing
    /// `c.4C>A` on NM_FOO.1 but requesting projection against NM_BAR.1)
    /// is rejected with a clear error rather than silently projecting
    /// through one transcript's exons then re-projecting through
    /// another's. For g. inputs there is no such input-transcript, so
    /// the caller's `transcript_id` is the only target.
    ///
    /// **3'-shift asymmetry note (re #334).** This routine receives a
    /// pre-normalized variant. The c./n./r. normalizer respects the
    /// HGVS exon-junction exception (does not shift across exons),
    /// while the g. normalizer does not. Consequently, the same
    /// biological variant fed as c. vs. g. near an exon junction can
    /// project to different `(g., c., p.)` tuples. This is intentional:
    /// the input axis carries the spec-correct semantic for *that*
    /// axis. Callers who want g.-axis semantics should pass the g.
    /// variant. The c.-input projection preserves the c.-axis canonical
    /// form by construction.
    /// (#328)
    /// Predict the protein consequence of a coding variant from its 1-based
    /// CDS position(s) and edit. Shared by the genome-pivot path
    /// (`project_single_inner`) and the direct bare-NM_ path
    /// (`project_coding_direct`) so both compute protein identically.
    ///
    /// Returns `Ok(None)` when no protein consequence applies — intronic, UTR,
    /// a non-coding transcript, or an edit shape with no in-frame/frameshift
    /// prediction. `Err` is reserved for genuine failures (unknown reference,
    /// etc.), never for "no consequence" (mirrors the existing contract).
    ///
    /// `cds_start`/`cds_end` are 1-based CDS positions; `cache_variant` keys
    /// the per-(transcript, parent) sequence and ref-translation caches.
    #[allow(clippy::too_many_arguments)]
    fn predict_protein_consequence(
        &self,
        transcript_id: &str,
        cdot_protein: Option<String>,
        is_coding: bool,
        is_intronic: bool,
        is_utr: bool,
        c_edit: &NaEdit,
        cds_start: &CdsPos,
        cds_end: &CdsPos,
        cache_variant: &HgvsVariant,
    ) -> Result<Option<HgvsVariant>, FerroError> {
        // An edit whose span reaches the translation initiation codon (CDS
        // 1–3) has an unpredictable protein consequence, reported as
        // `p.(Met1?)` (HGVS recommendations/protein/{substitution.md:51,
        // deletion.md:62}; #512). This holds even when the edit's 5' end lies
        // in the 5'UTR (CDS base ≤ 0) — e.g. a boundary-spanning insertion or
        // duplication like `c.-1_2dup` — so the coarse `is_utr` gate must not
        // drop it to "no protein predicted" (#504). Intronic offsets never
        // reach the initiation codon, so an offset disqualifies the edit here.
        let affects_init = !is_intronic
            && cds_start.offset.is_none()
            && cds_end.offset.is_none()
            && affects_initiation_codon(c_edit, cds_start.base, cds_end.base);
        if is_intronic || !is_coding || (is_utr && !affects_init) {
            return Ok(None);
        }
        // Resolve the protein accession: explicit cdot/reference value, else
        // RefSeq inference (NM_*→NP_*, XM_*→XP_*), else the transcript id
        // itself (the HGVS p. grammar requires a sequence identifier; see
        // #310). Substring stripping avoids mangling accessions whose suffix
        // contains `NM_`.
        let prot_acc: Option<String> = cdot_protein
            .or_else(|| {
                transcript_id
                    .strip_prefix("NM_")
                    .map(|rest| format!("NP_{rest}"))
            })
            .or_else(|| {
                transcript_id
                    .strip_prefix("XM_")
                    .map(|rest| format!("XP_{rest}"))
            })
            .or_else(|| Some(transcript_id.to_string()));
        let Some(prot_acc) = prot_acc else {
            return Ok(None);
        };
        // Issue #505: protein prediction below reads this transcript's CDS
        // bases directly and translates them. If the reference only carries a
        // *different* version of the transcript (a silent version-strip
        // substitution), those bases may pair with a different CDS/exon
        // structure, so the prediction would be wrong yet attributed to the
        // requested version. Decline to predict rather than emit a wrong
        // protein — this gates *all* downstream prediction, including the
        // initiation-codon short-circuit below, whose CDS positions are
        // likewise mapped against the requested-version transcript. The
        // predicate still permits exact-version cdot-genome synthesis — it
        // rejects only cross-version substitution.
        if !self.provider.has_transcript_version_exact(transcript_id) {
            return Ok(None);
        }
        // Issue #625: CDS-start sanity. Every protein prediction below trusts
        // that `Transcript.cds_start` points at the first base of the
        // initiation codon and translates the CDS bases directly from there.
        // When the cdot CDS annotation is inconsistent with the transcript
        // FASTA — cds_start lands mid-frame, not on an `ATG` (e.g.
        // NM_000425.3 cds_start→"AGG", NM_152263.2 cds_start→"GAT") — that
        // assumption is false: translation reads the wrong frame and would
        // fabricate a garbage consequence (premature internal stop, Xaa,
        // wrong anchor). Decline to predict rather than emit a bogus protein.
        //
        // The discriminator is the START codon only. Selenoproteins (e.g.
        // SELENON / NM_020451.2) legitimately carry an in-frame internal `TGA`
        // recoded as selenocysteine yet begin with `ATG`; keying on internal
        // stops would wrongly decline them, so we never inspect them here.
        //
        // This is a deliberate over-approximation: a few RefSeq transcripts use
        // experimentally-supported non-AUG initiation (CUG/GUG with a
        // `transl_except`), whose CDS is consistent yet non-`ATG`; we decline
        // those too rather than risk wrong-frame translation of the common
        // inconsistent-annotation case.
        //
        // A failure to read the CDS (missing sequence, degenerate coords) is
        // left to the per-edit paths below, which already treat
        // `ProteinSequenceUnavailable` as a non-fatal "no protein" — so we
        // only act on a CDS we could read whose start is not `ATG`. Read just
        // the start codon (not the whole 1–10 kbp CDS) since that is all the
        // guard inspects.
        let tx_for_cds_check =
            self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
        if let Ok(start_codon) = read_cds_start_codon(&tx_for_cds_check) {
            if !cds_has_valid_start(&start_codon) {
                log::trace!(
                    "declining protein prediction for {transcript_id}: reference CDS does \
                     not begin with an ATG start codon (cds_start is inconsistent with the \
                     transcript FASTA); first codon = {start_codon:?}",
                );
                return Ok(None);
            }
        }
        // An initiation-codon-affecting edit short-circuits to `p.(Met1?)`
        // here, before edit-type dispatch, so boundary-spanning edits whose
        // 5' end is in the 5'UTR (which the per-edit paths reject via their
        // `cds base > 0` guards) are still reported. The in-CDS start-codon
        // cases (e.g. `c.1A>G`) also funnel through here; the per-edit paths
        // keep their own `affects_initiation_codon` guard as defense in depth
        // for direct callers (#504, #512).
        if affects_init {
            let tx_for_codon =
                self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
            return Ok(Some(build_initiator_unknown(&prot_acc, &tx_for_codon)));
        }
        // Issue #332: route per-codon lookups through the variant-aware path
        // so an NG/NC-parented `cache_variant` picks the build-correct
        // chromosome; falls back to the bare provider lookup on
        // ReferenceNotFound. For a bare-NM_ `cache_variant` (no parent) this
        // is just the plain provider lookup.
        let mut protein = None;
        match c_edit {
            NaEdit::Substitution { .. } => {
                let tx_for_codon =
                    self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
                protein = Some(predict_substitution_protein(
                    &tx_for_codon,
                    cds_start.base,
                    c_edit,
                    &prot_acc,
                )?);
            }
            NaEdit::Deletion { .. }
            | NaEdit::Insertion { .. }
            | NaEdit::Duplication { .. }
            | NaEdit::Delins { .. }
            | NaEdit::Inversion { .. }
                // Only predict for concrete exonic CDS positions (intronic
                // offsets are already excluded by the is_intronic guard).
                if cds_start.offset.is_none()
                    && cds_end.offset.is_none()
                    && cds_start.base > 0
                    && cds_end.base > 0 =>
            {
                let tx_for_codon =
                    self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
                let ref_bundle =
                    self.cached_ref_translation(cache_variant, transcript_id, &tx_for_codon)?;
                match predict_indel_protein(
                    &tx_for_codon,
                    &ref_bundle,
                    cds_start.base,
                    cds_end.base,
                    c_edit,
                    &prot_acc,
                ) {
                    Ok(pv) => protein = Some(pv),
                    // Non-fatal: unsupported edits or missing sequence → None.
                    Err(FerroError::UnsupportedProjection { .. })
                    | Err(FerroError::ProteinSequenceUnavailable { .. }) => {}
                    Err(other) => return Err(other),
                }
            }
            _ => {}
        }
        Ok(protein)
    }

    /// Direct c.→p. projection for a bare transcript coding input (no
    /// `genomic_context` parent). The c.→g.→CDS roundtrip cannot run without a
    /// genome alignment, but protein prediction only needs the transcript's
    /// CDS sequence and the 1-based CDS position — which an exonic c. variant
    /// already provides. The resulting projection has `genomic = None` (no
    /// genomic representation is available for a bare-NM_ input) (#498).
    fn project_coding_direct(
        &self,
        cds: &CdsVariant,
        normalized: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Mirror the genome path's transcript_id-mismatch guard.
        let input_tx = cds.accession.transcript_accession();
        if input_tx != transcript_id {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "transcript_id mismatch: input is on {} but projection requested \
                     against {}; transcript-coordinate inputs must be projected against \
                     their own transcript",
                    input_tx, transcript_id,
                ),
            });
        }

        // Coding/protein/gene metadata: prefer cdot, fall back to the
        // sequence provider's transcript record.
        let (is_coding, cdot_protein, gene_symbol) =
            match self.projector.mapper().cdot().get_transcript(transcript_id) {
                Some(t) => (
                    t.cds_start.is_some(),
                    t.protein.clone(),
                    t.gene_name.clone(),
                ),
                None => {
                    let tx = self.cached_get_transcript_for_variant(normalized, transcript_id)?;
                    (
                        tx.cds_start.is_some(),
                        tx.protein_id.clone(),
                        tx.gene_symbol.clone(),
                    )
                }
            };

        let edit = cds.loc_edit.edit.inner().cloned().ok_or_else(|| {
            FerroError::UnsupportedProjection {
                reason: "coding variant has no concrete edit".to_string(),
            }
        })?;
        // Resolve the c. boundaries through the shared resolver so a bare `c.?`
        // (`Single(Mu::Unknown)`) surfaces as `UnsupportedProjection` — matching
        // the genome-pivot `project_to_genomic` contract — while a genuine range
        // boundary still surfaces as `InvalidCoordinates` (#508 review).
        let cds_start = resolve_uncertain_boundary(&cds.loc_edit.location.start, "c.", "start")?;
        let cds_end = resolve_uncertain_boundary(&cds.loc_edit.location.end, "c.", "end")?;

        // `resolve_uncertain_boundary` only rejects the parsed `?` sentinel
        // (`Single(Mu::Unknown)`). A programmatically-built `CdsPos::unknown`
        // carries its `base == 0` sentinel inside a `Mu::Certain`, so it slips
        // through above; without this guard `base <= 0` would misclassify it as
        // 5'UTR and return `Ok`. Reject it here to match the
        // `project_to_genomic` contract (#508 review).
        if cds_start.is_unknown()
            || cds_start.is_special()
            || cds_end.is_unknown()
            || cds_end.is_special()
        {
            return Err(FerroError::UnsupportedProjection {
                reason: "cannot resolve `?` or special (pter/qter/cen) position sentinels on \
                         direct c.→p. projection"
                    .to_string(),
            });
        }

        // Flags derived straight from the input c. position (no genome):
        // intronic = any offset; UTR = 5'UTR (base ≤ 0) or 3'UTR (`*`).
        let is_intronic = cds_start.offset.is_some() || cds_end.offset.is_some();
        let is_utr = !is_intronic
            && (cds_start.base <= 0 || cds_end.base <= 0 || cds_start.utr3 || cds_end.utr3);

        let protein = self.predict_protein_consequence(
            transcript_id,
            cdot_protein,
            is_coding,
            is_intronic,
            is_utr,
            &edit,
            &cds_start,
            &cds_end,
            normalized,
        )?;

        let frameshift = is_frameshift(normalized);

        // c↔n axis: derive the n. (transcript-relative) form from the resolved
        // CDS positions via the exon/CIGAR-aware mapper. Carry the *input's*
        // gene symbol (`cds.gene_symbol`) so the n. form matches the `coding`
        // form (which is the raw input here) rather than disagreeing on the
        // gene tag. Best-effort — a derivation edge (e.g. legacy c.0) is logged
        // and leaves the axis None rather than failing the whole projection.
        let noncoding = match self.cached_get_transcript_for_variant(normalized, transcript_id) {
            Ok(tx) => match crate::project::transcript_axis::noncoding_from_coding(
                &cds_start,
                &cds_end,
                &tx,
                &edit,
                transcript_id,
                cds.gene_symbol.clone(),
            ) {
                Ok(v) => Some(v),
                Err(e) => {
                    log::trace!("c→n derivation failed for {transcript_id}: {e}");
                    None
                }
            },
            Err(_) => None,
        };

        Ok(VariantProjection {
            genomic: None,
            coding: Some(normalized.clone()),
            noncoding,
            protein,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
        })
    }

    fn project_single_inner(
        &self,
        normalized: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Direct c.→p. path: a bare coding input (no genomic_context parent)
        // has no genome alignment to roundtrip through, but protein can be
        // predicted straight from the CDS (#498).
        if let HgvsVariant::Cds(c) = normalized {
            if c.accession.genomic_context.is_none() {
                return self.project_coding_direct(c, normalized, transcript_id);
            }
        }
        // Track which g. variant feeds the downstream pipeline AND will
        // be reported in `VariantProjection.genomic`. For g. input the
        // two are the same; for c./n./r. input we project once and use
        // the projected form for both — the `.genomic` field is the
        // canonical g. representation of the variant, not the input
        // axis.
        let projected_genome: HgvsVariant = match normalized {
            HgvsVariant::Genome(_) => normalized.clone(),
            HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_) => {
                // Reject transcript_id mismatch up front: projecting a
                // c.-on-NM_FOO through NM_BAR's exons and back through
                // NM_FOO's would silently produce nonsensical results.
                if let Some(acc) = normalized.accession() {
                    let input_tx = acc.transcript_accession();
                    if input_tx != transcript_id {
                        return Err(FerroError::UnsupportedProjection {
                            reason: format!(
                                "transcript_id mismatch: input is on {} but projection \
                                 requested against {}; transcript-coordinate inputs must \
                                 be projected against their own transcript",
                                input_tx, transcript_id,
                            ),
                        });
                    }
                }
                let g = self.project_to_genomic(normalized)?;
                match &g {
                    HgvsVariant::Genome(_) => g,
                    _ => {
                        // Defensive: project_to_genomic is documented to
                        // return Genome for these axes. If a future
                        // refactor changes that contract, surface a
                        // clear error rather than silently misprojecting.
                        return Err(FerroError::UnsupportedProjection {
                            reason: format!(
                                "project_to_genomic returned a non-Genome variant for {} input",
                                normalized.variant_type()
                            ),
                        });
                    }
                }
            }
            _ => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "VariantProjector does not accept {} inputs",
                        normalized.variant_type()
                    ),
                });
            }
        };
        let genome_variant = match &projected_genome {
            HgvsVariant::Genome(g) => g.clone(),
            _ => unreachable!("projected_genome is always Genome by construction above"),
        };

        let edit = genome_variant
            .loc_edit
            .edit
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::UnsupportedProjection {
                reason: "g. variant has no concrete edit".to_string(),
            })?;

        // Look up the transcript in the cdot mapper. When the projected
        // g. accession is NC_*, infer the genome build from its version
        // (NC_*.10 → GRCh37, NC_*.11 → GRCh38, etc.) so multi-build cdot
        // loads return the alignment matching the input chromosome rather
        // than silently using the primary build's view (issue #389).
        let build_hint =
            crate::liftover::aliases::infer_genome_build_from_accession(&genome_variant.accession);
        let cdot_tx = match build_hint {
            Some(b) => self
                .projector
                .mapper()
                .cdot()
                .get_transcript_on_build(transcript_id, b),
            None => self.projector.mapper().cdot().get_transcript(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;
        let strand = cdot_tx.strand;
        let gene_symbol = cdot_tx.gene_name.clone();
        let cdot_protein = cdot_tx.protein.clone();
        let is_coding = cdot_tx.cds_start.is_some();

        // Extract start and end genomic positions from the variant interval.
        let g_start = genome_variant
            .loc_edit
            .location
            .start
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "genomic interval start is unknown".to_string(),
            })?;
        let g_end = genome_variant
            .loc_edit
            .location
            .end
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "genomic interval end is unknown".to_string(),
            })?;

        let mapper = self.projector.mapper();

        // Genomic extent of the transcript — used to short-circuit the
        // genome-to-CDS mapping when the variant is outside the transcript's
        // span. The cached primary-build side-table is consulted via
        // `transcript_genome_span_on_build`; for alt builds the span is
        // computed on demand against the alt-build exon table (#389).
        let (tx_genome_start, tx_genome_end) = match build_hint {
            Some(b) => mapper
                .cdot()
                .transcript_genome_span_on_build(transcript_id, b),
            None => mapper.cdot().transcript_genome_span(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;

        // Helper: map one GenomePos → CdsPos, converting out-of-range errors.
        // `normalized.to_string()` is only consumed by the error message and
        // the happy path doesn't touch it; building it eagerly was ~2% of
        // SNP CPU per `fmt::write` self-time, so defer it to the moment we
        // actually construct a `TranscriptNotOverlapping` error.
        let map_position = |gp: &GenomePos| -> Result<(CdsPos, MappingInfo), FerroError> {
            if gp.base < tx_genome_start || gp.base >= tx_genome_end {
                return Err(FerroError::TranscriptNotOverlapping {
                    variant: normalized.to_string(),
                    transcript_id: transcript_id.to_string(),
                });
            }
            match mapper.genome_to_cds_on_build(transcript_id, gp, build_hint) {
                Ok(res) => Ok((res.variant, res.info)),
                Err(FerroError::InvalidCoordinates { .. })
                | Err(FerroError::ConversionError { .. }) => {
                    Err(FerroError::TranscriptNotOverlapping {
                        variant: normalized.to_string(),
                        transcript_id: transcript_id.to_string(),
                    })
                }
                Err(other) => Err(other),
            }
        };

        let (cds_start_raw, info_start) = map_position(&g_start)?;
        let (cds_end_raw, info_end) = map_position(&g_end)?;

        // On minus strand the start and end of the c. interval are swapped.
        let (cds_start, cds_end) = match strand {
            crate::reference::Strand::Plus => (cds_start_raw, cds_end_raw),
            crate::reference::Strand::Minus => (cds_end_raw, cds_start_raw),
            crate::reference::Strand::Unknown => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "transcript {} has unknown strand; cannot project g. → c./n.",
                        transcript_id
                    ),
                });
            }
        };

        // Transform the edit for the transcript strand.
        let c_edit = transform_edit_for_strand(&edit, strand);

        // Build the c./n. HGVS variant returned in `VariantProjection.coding`.
        // Kept bare-accession (no NC_* wrapper) for caller-visible rendering;
        // build-awareness for the parent-aware caches goes through
        // `cache_variant` below.
        let coding = if is_coding {
            let interval = CdsInterval::new(cds_start, cds_end);
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession(transcript_id),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, c_edit.clone()),
            })
        } else {
            let tx_start = TxPos::new(cds_start.base);
            let tx_end = TxPos::new(cds_end.base);
            HgvsVariant::Tx(TxVariant {
                accession: parse_accession(transcript_id),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(TxInterval::new(tx_start, tx_end), c_edit.clone()),
            })
        };

        // Parallel "cache key" variant that stamps the projected g.
        // accession (NC_*.10 / NC_*.11) as `genomic_context`. Used only
        // for `cached_get_transcript_for_variant` and
        // `cached_ref_translation` so that:
        //   1. `parent_key_for` reads the build-bearing NC accession,
        //      partitioning the (transcript_id, parent) caches by
        //      build — NC_*.10 vs NC_*.11 inputs no longer collide
        //      under the same key.
        //   2. `MultiFastaProvider::get_transcript_for_variant` sees the
        //      NC parent and probes the matching cdot build first via
        //      `get_transcript_on_build`, so the cached `Transcript` and
        //      `RefProteinBundle` are built against the correct
        //      alignment.
        // The returned `coding` field is unaffected (#389 follow-up).
        let cache_variant = if is_coding {
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession(transcript_id)
                    .with_genomic_context(genome_variant.accession.clone()),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(CdsInterval::new(cds_start, cds_end), c_edit.clone()),
            })
        } else {
            HgvsVariant::Tx(TxVariant {
                accession: parse_accession(transcript_id)
                    .with_genomic_context(genome_variant.accession.clone()),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(
                    TxInterval::new(TxPos::new(cds_start.base), TxPos::new(cds_end.base)),
                    c_edit.clone(),
                ),
            })
        };

        // Derive per-position flags from MappingInfo.
        let is_intronic = info_start.is_intronic || info_end.is_intronic;
        let is_utr = !is_intronic
            && (info_start.in_5utr || info_start.in_3utr || info_end.in_5utr || info_end.in_3utr);

        // Predict protein consequence for CDS variants. Shared with the
        // direct bare-NM_ c.→p. path (`project_coding_direct`) so both entry
        // points use identical CDS→protein logic (#498).
        let protein = self.predict_protein_consequence(
            transcript_id,
            cdot_protein,
            is_coding,
            is_intronic,
            is_utr,
            &c_edit,
            &cds_start,
            &cds_end,
            &cache_variant,
        )?;

        let frameshift = is_frameshift(&coding);

        // c↔n axis. For a coding transcript derive the n. form genome-free from
        // the resolved CDS positions (same edit, reframed coordinates);
        // best-effort — never fail the whole projection on an n.-derivation
        // edge. For a non-coding transcript `coding` already IS the n. (Tx)
        // form, so reuse it (forward-compatible: the non-coding full-projection
        // path is currently pinned unsupported — see
        // `tests/issue_328_projector_accepts_cnr.rs` — so this arm is not yet
        // reached, but carries the correct semantics for when it is).
        let noncoding = if is_coding {
            match self.cached_get_transcript_for_variant(&cache_variant, transcript_id) {
                Ok(tx) => match crate::project::transcript_axis::noncoding_from_coding(
                    &cds_start,
                    &cds_end,
                    &tx,
                    &c_edit,
                    transcript_id,
                    gene_symbol.clone(),
                ) {
                    Ok(v) => Some(v),
                    Err(e) => {
                        log::trace!("c→n derivation failed for {transcript_id}: {e}");
                        None
                    }
                },
                Err(_) => None,
            }
        } else {
            Some(coding.clone())
        };

        // `.genomic` is always the canonical g. representation. For g.
        // input that's the (normalized) input itself; for c./n./r.
        // input that's the variant produced by `project_to_genomic`.
        // Using `projected_genome` for both cases keeps `.genomic`
        // axis-correct regardless of how the caller entered.
        Ok(VariantProjection {
            genomic: Some(projected_genome),
            coding: Some(coding),
            noncoding,
            protein,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
        })
    }
}

/// Compute a representative 1-based genomic position for a non-coding
/// (n.) or RNA (r.) transcript coordinate, suitable for seeding
/// `Projector::project`'s stab query.
///
/// - For "simple" exonic positions (`offset.is_none() && !utr3 &&
///   base >= 1`), the exact `tx_to_genome` mapping is returned.
/// - For intronic (`offset.is_some()`), downstream/3'UTR (`utr3`), or
///   non-positive `base`, the conversion is best-effort: the
///   transcript's first-exon genome_start is returned. The stab query
///   at that position still lands inside the transcript, and the
///   per-transcript projection downstream re-derives the precise
///   coordinate. Falling back here keeps the fan-out path lenient for
///   inputs the precise tx-to-genome path can't resolve (which is also
///   what `project_to_genomic` itself rejects elsewhere for intronic
///   non-coding offsets).
fn nr_representative_genome_pos(
    cdot_tx: &CdotTranscript,
    base: i64,
    offset: Option<i64>,
    utr3: bool,
) -> Result<u64, FerroError> {
    if offset.is_some() || utr3 || base < 1 {
        return cdot_tx
            .exons
            .first()
            .map(|e| e[0])
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "transcript has no exons; cannot derive stab position".to_string(),
            });
    }
    let tx_pos_0based = (base - 1) as u64;
    cdot_tx
        .tx_to_genome(tx_pos_0based)
        .ok_or_else(|| FerroError::InvalidCoordinates {
            msg: format!(
                "tx position {} not in any exon of {} (cannot derive stab position)",
                base, cdot_tx.contig
            ),
        })
}

/// Resolve an [`UncertainBoundary<T>`] to an owned `T`, distinguishing the two
/// `None` returns of [`UncertainBoundary::inner`]:
///   - `Single(Mu::Unknown)` (parsed `?` sentinel) → [`FerroError::UnsupportedProjection`]
///   - `Range { .. }`        (parsed `(a_b)` range) → [`FerroError::InvalidCoordinates`]
///
/// Without this split, a parsed `c.?` slips through as `InvalidCoordinates` and
/// contradicts the documented projection contract. Shared by both projection
/// entry points — `project_to_genomic` (genome-pivot path) and
/// `project_coding_direct` (direct bare-NM_ path) — so they classify `?` vs.
/// ranges identically.
fn resolve_uncertain_boundary<T: Copy>(
    boundary: &crate::hgvs::interval::UncertainBoundary<T>,
    coord_label: &str,
    end_label: &str,
) -> Result<T, FerroError> {
    use crate::hgvs::interval::UncertainBoundary;
    use crate::hgvs::uncertainty::Mu;
    match boundary {
        UncertainBoundary::Single(Mu::Certain(t) | Mu::Uncertain(t)) => Ok(*t),
        UncertainBoundary::Single(Mu::Unknown) => Err(FerroError::UnsupportedProjection {
            reason: format!("cannot resolve unknown position (`?`) at {coord_label} {end_label}"),
        }),
        UncertainBoundary::Range { .. } => Err(FerroError::InvalidCoordinates {
            msg: format!("{coord_label} interval {end_label} is a range"),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::{CdotMapper, CdotTranscript};
    use crate::data::projection::Projector;
    use crate::hgvs::variant::{Accession, AllelePhase, AlleleVariant};
    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
    use crate::reference::Strand as ProvStrand;
    use std::sync::OnceLock;

    fn make_test_provider_and_projector() -> (Projector, MockProvider) {
        // A single coding transcript on chr1, plus strand.
        //
        // Genomic layout (cdot 0-based coords):
        //   Exon: genome [1000, 1009), tx [0, 9), so 9 bases total.
        //   cds_start = 0 (0-based, no 5'UTR), cds_end = 9.
        //
        // Sequence: "ATGCGCTAA" = Met-Arg-Stop (3 codons).
        //
        // Coordinate mapping (cdot genome 0-based → HGVS c. 1-based):
        //   g.1000 → tx_pos=0 → c.1 (A = first base of Met)
        //   g.1001 → tx_pos=1 → c.2 (T)
        //   g.1002 → tx_pos=2 → c.3 (G)
        //   g.1003 → tx_pos=3 → c.4 (C ← ref base for test substitution)
        //   g.1004 → tx_pos=4 → c.5 (G)
        //   ...
        //
        // c.4C>A: codon 2 CGC (Arg) → AGC (Ser), missense.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_TEST.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // [genome_start(0-based), genome_end(0-based excl), tx_start(1-based), tx_end(1-based)]
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0), // 0-based: CDS starts at tx_pos 0 (no 5'UTR)
                cds_end: Some(9),   // 0-based exclusive: CDS ends at tx_pos 9
                gene_id: None,
                protein: Some("NP_TEST.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        // Transcript: sequence "ATGCGCTAA", cds_start=1 (1-based, first base).
        provider.add_transcript(Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCTAA".to_string()),
            cds_start: Some(1), // 1-based inclusive per Transcript convention
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Genomic sequence: 999 N's + "ATGCGCTAA" + 100 N's.
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    fn make_minus_strand_provider_and_projector() -> (Projector, MockProvider) {
        // Same 9bp CDS on chr1, but transcript is on the minus strand.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_TEST_MINUS.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Minus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_TEST_MINUS.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_TEST_MINUS.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Minus,
            sequence: Some("ATGCGCTAA".to_string()), // CDS as the transcript reads it
            cds_start: Some(1),                      // 1-based inclusive per Transcript convention
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}TTAGCGCAT{}", prefix, suffix));
        (projector, provider)
    }

    /// Build a two-transcript setup for project_all tests.
    ///
    /// NM_TX1.1: chr1 [1000,1009), plus strand, 9bp CDS "ATGCGCTAA"
    /// NM_TX2.1: chr1 [1000,1009), plus strand, 9bp CDS "ATGCGCTAA" (same region)
    ///           NM_TX2.1 is registered as MANE Select so it sorts first.
    fn make_two_transcript_setup() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_TX1.1".to_string(),
            CdotTranscript {
                gene_name: Some("GENE1".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_TX1.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        cdot.add_transcript(
            "NM_TX2.1".to_string(),
            CdotTranscript {
                gene_name: Some("GENE1".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_TX2.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot)
            // Make NM_TX2.1 the MANE Select so it sorts first.
            .with_mane(vec!["NM_TX2.1".to_string()], vec![]);

        let mut provider = MockProvider::new();
        for id in ["NM_TX1.1", "NM_TX2.1"] {
            provider.add_transcript(Transcript {
                id: id.to_string(),
                gene_symbol: Some("GENE1".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAA".to_string()),
                cds_start: Some(1),
                cds_end: Some(9),
                exons: vec![Exon::new(1, 1, 9)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1008),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
        }
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    // -------------------------------------------------------------------------
    // Cache identity tests
    // -------------------------------------------------------------------------

    /// Build a minimal coding HgvsVariant pointing at `transcript_id` for use
    /// in `cached_ref_translation` tests. The variant's `accession` carries
    /// `genomic_context = ctx`, so `parent_key_for` returns the same string
    /// the production code would derive at projection time.
    #[cfg(test)]
    fn make_coding_variant(transcript_id: &str, ctx: Option<Accession>) -> HgvsVariant {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};

        let mut accession = parse_accession(transcript_id);
        if let Some(parent) = ctx {
            accession = accession.with_genomic_context(parent);
        }
        HgvsVariant::Cds(CdsVariant {
            accession,
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        })
    }

    /// Two consecutive `cached_ref_translation` calls for the same
    /// `(transcript_id, parent_accession)` pair must return the same
    /// `Arc<RefProteinBundle>` (pointer-equal). Catches any future
    /// regression that accidentally rebuilds the bundle per call.
    #[test]
    fn cached_ref_translation_returns_same_arc() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let tx = vp
            .cached_get_transcript("NM_TEST.1")
            .expect("transcript fetch");
        let variant = make_coding_variant("NM_TEST.1", None);
        let a = vp
            .cached_ref_translation(&variant, "NM_TEST.1", &tx)
            .expect("first translation");
        let b = vp
            .cached_ref_translation(&variant, "NM_TEST.1", &tx)
            .expect("second translation");
        assert!(
            Arc::ptr_eq(&a, &b),
            "ref-protein cache must return the same Arc for repeated lookups"
        );
        // ATGCGCTAA = Met-Arg-Ter; translate_full_cds stops before the Ter.
        assert_eq!(a.ref_protein.len(), 2);
        assert_eq!(a.ref_protein_with_stop.len(), 3);
    }

    /// Two `cached_ref_translation` lookups for the same `transcript_id`
    /// but different parent `genomic_context` accessions must produce
    /// distinct cache entries (different `Arc`s). Without parent-aware
    /// keying the second call would reuse the bundle built for the first
    /// parent, which can poison indel protein predictions once the
    /// underlying transcript resolves to a different cdot build under
    /// each parent (issue #332).
    #[test]
    fn cached_ref_translation_is_parent_aware() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let tx = vp
            .cached_get_transcript("NM_TEST.1")
            .expect("transcript fetch");

        let v_none = make_coding_variant("NM_TEST.1", None);
        let v_ng1 = make_coding_variant("NM_TEST.1", Some(Accession::new("NG", "TEST", Some(1))));
        let v_ng2 = make_coding_variant("NM_TEST.1", Some(Accession::new("NG", "TEST", Some(2))));

        let b_none = vp
            .cached_ref_translation(&v_none, "NM_TEST.1", &tx)
            .expect("no-parent bundle");
        let b_ng1 = vp
            .cached_ref_translation(&v_ng1, "NM_TEST.1", &tx)
            .expect("NG_TEST.1 bundle");
        let b_ng2 = vp
            .cached_ref_translation(&v_ng2, "NM_TEST.1", &tx)
            .expect("NG_TEST.2 bundle");

        // Each parent identity must occupy its own cache slot.
        assert!(
            !Arc::ptr_eq(&b_none, &b_ng1),
            "no-parent and NG-parented entries must not collide"
        );
        assert!(
            !Arc::ptr_eq(&b_ng1, &b_ng2),
            "two distinct NG parents must produce distinct cache entries"
        );

        // Re-querying with the same parent must still hit the cache.
        let b_ng1_again = vp
            .cached_ref_translation(&v_ng1, "NM_TEST.1", &tx)
            .expect("NG_TEST.1 bundle (repeat)");
        assert!(
            Arc::ptr_eq(&b_ng1, &b_ng1_again),
            "repeated lookup with the same parent must hit the cache"
        );
    }

    // -------------------------------------------------------------------------
    // Existing tests (preserved)
    // -------------------------------------------------------------------------

    #[test]
    fn project_substitution_minus_strand_revcomps_ref_alt() {
        let (projector, provider) = make_minus_strand_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1005G>A", "NM_TEST_MINUS.1")
            .expect("minus-strand projection should succeed");
        let c = result
            .coding
            .as_ref()
            .expect("c. should be present")
            .to_string();
        assert!(
            c.contains(":c.4C>T"),
            "expected revcomp c.4C>T for minus strand, got: {}",
            c
        );
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present")
            .to_string();
        assert_eq!(p, "NP_TEST_MINUS.1:p.(Arg2Cys)");
        assert!(!result.is_frameshift);
        assert!(!result.is_intronic);
        assert!(!result.is_utr);
    }

    fn make_intronic_test_data() -> (Projector, MockProvider) {
        // Two-exon coding transcript on chr1, plus strand.
        //
        //   Exon 1: genome [1000, 1010), tx [0, 10)  — 10 bp
        //   Intron: genome [1010, 2000)               — 990 bp
        //   Exon 2: genome [2000, 2010), tx [10, 20) — 10 bp (last 2 are pad)
        //
        // CDS: tx [0, 18) → 18 bp → 6 codons: ATG-CGC-AAA-GGG-TAA-CCC
        //   (Met-Arg-Lys-Gly-Stop-Pro)
        //
        // Test position g.1015 is in the intron, 6 bases after the end of exon 1
        // (exon 1 ends at genome 1010, exclusive; last exonic base is g.1009;
        //  distance = 1015-1010+1 = 6).  Maps to c.9+6 or c.10+6 depending on
        // whether tx_end is inclusive or exclusive in the boundary calculation.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_INTR.1".to_string(),
            CdotTranscript {
                gene_name: Some("INTRGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1010, 0, 10], [2000, 2010, 10, 20]],
                cds_start: Some(0),
                cds_end: Some(18),
                gene_id: None,
                protein: Some("NP_INTR.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_INTR.1".to_string(),
            gene_symbol: Some("INTRGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCAAAGGGTAACCC".to_string()), // 18 bp
            cds_start: Some(1),
            cds_end: Some(18),
            exons: vec![Exon::new(1, 1, 10), Exon::new(2, 11, 20)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2009),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let exon1 = "ATGCGCAAAG"; // 10 bp (genomic positions 1000..1009)
        let intron = "N".repeat(990); // genomic positions 1010..1999 (990 bp)
        let exon2 = "GGTAACCCNN"; // 10 bp pad (genomic positions 2000..2009)
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{prefix}{exon1}{intron}{exon2}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_intronic_substitution_no_protein() {
        let (projector, provider) = make_intronic_test_data();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1015A>G", "NM_INTR.1")
            .expect("intronic substitution should project to c. with offset");
        assert!(result.is_intronic, "expected is_intronic=true");
        assert!(result.protein.is_none(), "no p. for intronic substitutions");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(
            c.contains('+'),
            "expected intronic offset notation (e.g. c.9+6 or c.10+6), got: {}",
            c
        );
    }

    #[test]
    fn project_no_overlap_returns_transcript_not_overlapping() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let err = vp
            .project("NC_000001.11:g.5000A>G", "NM_TEST.1")
            .expect_err("should fail to project outside the transcript");
        assert!(
            matches!(err, FerroError::TranscriptNotOverlapping { .. }),
            "expected TranscriptNotOverlapping, got: {:?}",
            err
        );
    }

    #[test]
    fn project_substitution_plus_strand_missense() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        let result = vp
            .project("chr1:g.1003C>A", "NM_TEST.1")
            .expect("projection should succeed");

        assert_eq!(result.transcript_id, "NM_TEST.1");
        assert_eq!(result.gene_symbol.as_deref(), Some("TESTGENE"));

        let c = result
            .coding
            .as_ref()
            .expect("c. should be present")
            .to_string();
        assert!(c.contains(":c.4C>A"), "expected ':c.4C>A' in '{}' ", c);

        let p = result
            .protein
            .as_ref()
            .expect("p. should be present")
            .to_string();
        assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");

        assert!(!result.is_frameshift);
        assert!(!result.is_intronic);
        assert!(!result.is_utr);
    }

    // -------------------------------------------------------------------------
    // CDS-start sanity (issue #625)
    //
    // When a transcript's cdot CDS coordinates are inconsistent with its FASTA
    // (cds_start does not land on an ATG start codon), the reference CDS reads
    // the wrong frame. ferro must DECLINE protein prediction (Ok(None)) rather
    // than emit a fabricated consequence translated from the wrong frame.
    // -------------------------------------------------------------------------

    /// Build a single-transcript projector whose CDS sequence and length are
    /// chosen by the caller. The whole transcript is one plus-strand exon at
    /// genome `[1000, 1000+len)`, with `cds_start = 1` and `cds_end = len`
    /// (1-based, full-length CDS, no UTR). The genomic sequence mirrors the
    /// transcript so a `g.(1000+k)` SNV maps to `c.(k+1)`.
    ///
    /// `seq` is the transcript/CDS sequence; its first codon is what the
    /// CDS-start sanity guard inspects.
    #[cfg(test)]
    fn make_custom_cds_provider_and_projector(seq: &str) -> (Projector, MockProvider) {
        let len_u = seq.len() as u64;
        let genome_end_excl_u = 1000_u64 + len_u; // 0-based exclusive

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_CUSTOM.1".to_string(),
            CdotTranscript {
                gene_name: Some("CUSTOMGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end]
                exons: vec![[1000, genome_end_excl_u, 0, len_u]],
                cds_start: Some(0),   // 0-based: CDS starts at tx_pos 0 (no 5'UTR)
                cds_end: Some(len_u), // 0-based exclusive
                gene_id: None,
                protein: Some("NP_CUSTOM.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_CUSTOM.1".to_string(),
            gene_symbol: Some("CUSTOMGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some(seq.to_string()),
            cds_start: Some(1), // 1-based inclusive per Transcript convention
            cds_end: Some(len_u),
            exons: vec![Exon::new(1, 1, len_u)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(genome_end_excl_u - 1),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(1000);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{prefix}{seq}{suffix}"));
        (projector, provider)
    }

    /// (a) RED-FIRST: a transcript whose `cds_start` lands on a non-ATG codon
    /// ("AGG", not a start) is FASTA/CDS-inconsistent. A coding SNV must
    /// DECLINE protein prediction (no fabricated consequence from the wrong
    /// frame), even though c. and n. axes are still emitted.
    #[test]
    fn project_non_atg_start_declines_protein() {
        // CDS "AGGGGCGCTAA": first codon AGG (Arg) is NOT a start codon.
        // g.1003 → c.4 (4th base 'G'), a plainly coding position.
        let (projector, provider) = make_custom_cds_provider_and_projector("AGGGGCGCTAA");
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003G>A", "NM_CUSTOM.1")
            .expect("projection should still succeed at the c./n. level");

        // c. axis must still be present (the inconsistency is CDS-only).
        assert!(
            result.coding.is_some(),
            "c. axis should still be emitted for a coding SNV"
        );
        // Protein must be DECLINED: the CDS does not start with ATG, so any
        // translation would be from the wrong frame.
        assert!(
            result.protein.is_none(),
            "expected NO protein for a non-ATG-start (inconsistent) CDS, got: {:?}",
            result.protein
        );
    }

    /// (b) CONTROL: a normal ATG-start CDS predicts protein as before. This
    /// must stay green after the guard is added.
    #[test]
    fn project_atg_start_predicts_protein_control() {
        // CDS "ATGCGCTAA": Met-Arg-Ter, a valid ATG start.
        let (projector, provider) = make_custom_cds_provider_and_projector("ATGCGCTAA");
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "NM_CUSTOM.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present for a valid ATG-start CDS")
            .to_string();
        assert_eq!(p, "NP_CUSTOM.1:p.(Arg2Ser)");
    }

    /// (c) SELENOPROTEIN CONTROL: a valid ATG-start CDS that contains an
    /// in-frame internal TGA (recoded as selenocysteine in real selenoproteins
    /// such as SELENON / NM_020451.2) before the true stop. The guard MUST key
    /// on the START codon only — an internal stop is legitimate, so protein
    /// prediction must still succeed (NOT be declined).
    #[test]
    fn project_selenoprotein_internal_tga_predicts_protein() {
        // CDS "ATGAAATGAAAATAA":
        //   codon0 ATG (Met, start)  codon1 AAA (Lys)  codon2 TGA (internal,
        //   Sec in vivo / Ter for the naive translator)  codon3 AAA (Lys)
        //   codon4 TAA (true Ter).
        // A naive translator stops at the internal TGA, but the START is a
        // valid ATG, so the CDS-start guard must NOT decline.
        let (projector, provider) = make_custom_cds_provider_and_projector("ATGAAATGAAAATAA");
        let vp = VariantProjector::new(projector, provider);
        // g.1004 → c.5: 2nd base of codon1 (AAA, Lys) — squarely coding and
        // upstream of the internal TGA so prediction is well-defined.
        let result = vp
            .project("chr1:g.1004A>G", "NM_CUSTOM.1")
            .expect("projection should succeed");
        // Assert the concrete consequence, not just presence: codon1 AAA (Lys2)
        // → AGA (Arg2). A regression that declined selenoproteins — or emitted a
        // wrong-but-present protein — would slip past a bare is_some() check.
        let p = result
            .protein
            .as_ref()
            .expect(
                "selenoprotein-style internal TGA must NOT decline protein \
                 (guard keys on START codon, not internal stops)",
            )
            .to_string();
        assert_eq!(p, "NP_CUSTOM.1:p.(Lys2Arg)");
    }

    #[test]
    fn project_genome_pivot_coding_populates_noncoding_axis() {
        // c↔n axis: a g. SNV projected onto a coding transcript carries BOTH
        // the c. form and the n. form, in distinct coordinate frames.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "NM_TEST.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present");
        assert!(
            matches!(c, HgvsVariant::Cds(_)),
            "coding should be a c. (Cds) variant, got {c:?}"
        );
        let n = result
            .noncoding
            .as_ref()
            .expect("n. axis populated for coding tx");
        assert!(
            matches!(n, HgvsVariant::Tx(_)),
            "noncoding should be an n. (Tx) variant, got {n:?}"
        );
        // NM_TEST.1's CDS begins at transcript position 1, so c.4 → n.4.
        assert!(
            n.to_string().contains(":n.4C>A"),
            "expected n.4C>A, got {n}"
        );
    }

    #[test]
    fn project_single_base_deletion_is_frameshift() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1004del", "NM_TEST.1")
            .expect("deletion should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains(":c."), "expected c. variant, got: {}", c);
        assert!(c.contains("del"), "expected del notation in c., got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for CDS del")
            .to_string();
        assert!(p.contains("Arg2"), "expected Arg2 in p.: {}", p);
        assert!(p.contains("fs"), "expected fs in p.: {}", p);
        assert!(result.is_frameshift, "1-base del should be frameshift");
        assert!(!result.is_intronic);
    }

    #[test]
    fn project_three_base_deletion_in_frame() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1005del", "NM_TEST.1")
            .expect("3-base del should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("del"), "expected del notation, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for in-frame del")
            .to_string();
        assert!(p.contains("Arg2del"), "expected Arg2del in p.: {}", p);
        assert!(!result.is_frameshift, "3-base deletion is in-frame");
    }

    #[test]
    fn project_single_base_insertion_is_frameshift() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1004insA", "NM_TEST.1")
            .expect("insertion should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("ins"), "expected ins notation, got: {}", c);
        assert!(result.protein.is_some(), "p. expected for CDS insertion");
        assert!(result.is_frameshift, "1-base insertion is frameshift");
    }

    #[test]
    fn project_duplication() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1004dup", "NM_TEST.1")
            .expect("dup should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("dup"), "expected dup notation, got: {}", c);
        assert!(result.protein.is_some(), "p. expected for CDS dup");
        assert!(result.is_frameshift, "1-base dup is frameshift");
    }

    #[test]
    fn project_delins() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1004delinsAT", "NM_TEST.1")
            .expect("delins should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("delins"), "expected delins notation, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for delins")
            .to_string();
        assert!(p.contains("delins"), "expected delins in p.: {}", p);
        assert!(
            !result.is_frameshift,
            "delins of equal length is not a frameshift"
        );
    }

    #[test]
    fn project_inversion() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1005inv", "NM_TEST.1")
            .expect("inv should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("inv"), "expected inv notation, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for inv")
            .to_string();
        assert!(p.contains("Arg2"), "expected Arg2 in p.: {}", p);
        assert!(p.contains("delins"), "expected delins in p.: {}", p);
        assert!(p.contains("Ala"), "expected Ala in p.: {}", p);
        assert!(!result.is_frameshift, "inversion is not a frameshift");
    }

    fn make_ensembl_provider_and_projector() -> (Projector, MockProvider) {
        // Same CDS as make_test_provider_and_projector but the transcript ID is
        // an Ensembl accession (no NM_/XM_ prefix) and cdot.protein is absent.
        // The projector must not fabricate a bogus accession via substring
        // substitution; per #310 it falls back to using the transcript id as
        // the `p.` accession prefix.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "ENST00000000001.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "ENST00000000001.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    #[test]
    fn project_substitution_falls_back_to_transcript_id_as_protein_accession() {
        // Per #310: when no cdot.protein and no NM_/XM_ prefix is available
        // to infer NP_/XP_ from, the projector falls back to using the
        // transcript id directly as the `p.` accession prefix rather than
        // silently dropping the protein prediction. The HGVS protein
        // grammar requires a sequence_identifier, so the alternative would
        // be a non-conformant accession-less `p.` form.
        let (projector, provider) = make_ensembl_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "ENST00000000001.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains(":c.4C>A"), "expected c.4C>A, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present via transcript-id fallback")
            .to_string();
        assert_eq!(p, "ENST00000000001.1:p.(Arg2Ser)");
    }

    // -------------------------------------------------------------------------
    // project_normalized tests
    // -------------------------------------------------------------------------

    #[test]
    fn project_normalized_same_result_as_project_variant() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        let variant = crate::parse_hgvs("chr1:g.1003C>A").expect("parse should succeed");
        let via_project = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("project_variant should succeed");
        // Pre-normalize via the cached normalizer, then call project_normalized.
        let normalized = vp.normalizer.normalize(&variant).expect("normalize failed");
        let via_normalized = vp
            .project_normalized(&normalized, "NM_TEST.1")
            .expect("project_normalized should succeed");

        assert_eq!(
            via_project.coding.as_ref().map(|v| v.to_string()),
            via_normalized.coding.as_ref().map(|v| v.to_string()),
            "project_normalized should produce same c. as project_variant"
        );
        assert_eq!(
            via_project.protein.as_ref().map(|v| v.to_string()),
            via_normalized.protein.as_ref().map(|v| v.to_string()),
            "project_normalized should produce same p. as project_variant"
        );
    }

    // -------------------------------------------------------------------------
    // project_all / project_variant_all tests
    // -------------------------------------------------------------------------

    #[test]
    fn project_variant_all_returns_both_transcripts() {
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("chr1:g.1003C>A")
            .expect("project_all should succeed");

        assert_eq!(
            results.len(),
            2,
            "expected projections onto both transcripts, got: {}",
            results.len()
        );
    }

    #[test]
    fn project_all_mane_select_sorts_first() {
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("chr1:g.1003C>A")
            .expect("project_all should succeed");

        // NM_TX2.1 is MANE Select → must be first.
        assert_eq!(
            results[0].transcript_id, "NM_TX2.1",
            "MANE Select should be first"
        );
    }

    // -------------------------------------------------------------------------
    // c./n./r. fan-out via project_variant_all (closes #389 item 3)
    //
    // Pre-fix, `extract_contig_and_pos` matched only `HgvsVariant::Genome`
    // and returned `UnsupportedProjection { reason: "project_all
    // currently only accepts g. variants" }` for every other axis, which
    // defeated PR #379's widening of `project()` to c./n./r. inputs in
    // the fan-out entrypoints.
    // -------------------------------------------------------------------------

    /// `project_variant_all` on a c. input with an NG/NC parent must
    /// reach the fan-out path and produce a per-transcript projection
    /// for the matching transcript. Pre-#389 this errored out at the
    /// `extract_contig_and_pos` g.-only gate.
    #[test]
    fn project_variant_all_accepts_coding_input_with_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, HgvsVariant, LocEdit};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // NM_TX1.1 c.4C>A with an NC_000001.11 parent; per the
        // make_two_transcript_setup fixture this should project to
        // g.1003 (which is c.4 on the plus-strand 9-base CDS).
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let results = vp
            .project_variant_all(&HgvsVariant::Cds(cds))
            .expect("project_variant_all should accept c. input post-#389");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected at least the input transcript's projection, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// Same as above but for the n. axis: `project_variant_all` on a
    /// non-coding (or n.-axis) variant with an NG/NC parent must also
    /// reach the fan-out path.
    #[test]
    fn project_variant_all_accepts_noncoding_n_input_with_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{Accession, HgvsVariant, LocEdit, TxVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // n.4 maps to genome 1003 on NM_TX1.1 in this fixture (tx_pos 3
        // → exon.genome_start + 3 = 1003).
        let tx = TxVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let results = vp
            .project_variant_all(&HgvsVariant::Tx(tx))
            .expect("project_variant_all should accept n. input post-#389");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected at least the input transcript's projection, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// `project_variant_all` on an r. (RNA) input — parallel to n. but
    /// on the r. axis. Pre-#389 the gate rejected this too.
    #[test]
    fn project_variant_all_accepts_rna_input_with_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{Accession, HgvsVariant, LocEdit, RnaVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let rna = RnaVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let results = vp
            .project_variant_all(&HgvsVariant::Rna(rna))
            .expect("project_variant_all should accept r. input post-#389");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected at least the input transcript's projection, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// c. input *without* a `genomic_context` parent: the fan-out
    /// rejects the input up front with `UnsupportedProjection`
    /// (#389 follow-up).
    ///
    /// Pre-fix the stab query still seeded and each per-transcript
    /// projection silently failed inside `project_to_genomic`'s parent
    /// requirement; the fan-out loop's `log::trace!`-and-drop converted
    /// the user error into `Ok([])`, turning "you forgot the parent"
    /// into "no overlaps." Post-fix the error is surfaced.
    #[test]
    fn project_variant_all_rejects_coding_input_without_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{CdsVariant, HgvsVariant, LocEdit};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // No `.with_genomic_context(...)` — bare NM accession on a c. input.
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1"),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let err = vp
            .project_variant_all(&HgvsVariant::Cds(cds))
            .expect_err("c. input without parent must error, not silently return empty");
        match err {
            FerroError::UnsupportedProjection { reason } => {
                assert!(
                    reason.contains("project_all requires a parent reference")
                        && reason.contains("c."),
                    "unexpected error message: {reason}"
                );
            }
            other => panic!("expected UnsupportedProjection, got {other:?}"),
        }
    }

    /// Same as above but for n. — the Tx-axis arm of the
    /// `require_parent_for_fanout` gate (#389 follow-up).
    #[test]
    fn project_variant_all_rejects_noncoding_n_input_without_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{HgvsVariant, LocEdit, TxVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let tx = TxVariant {
            accession: parse_accession("NM_TX1.1"),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let err = vp
            .project_variant_all(&HgvsVariant::Tx(tx))
            .expect_err("n. input without parent must error");
        match err {
            FerroError::UnsupportedProjection { reason } => {
                assert!(reason.contains("n."), "unexpected error message: {reason}");
            }
            other => panic!("expected UnsupportedProjection, got {other:?}"),
        }
    }

    /// Same as above but for r. — the Rna-axis arm of the gate.
    #[test]
    fn project_variant_all_rejects_rna_input_without_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{HgvsVariant, LocEdit, RnaVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let rna = RnaVariant {
            accession: parse_accession("NM_TX1.1"),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let err = vp
            .project_variant_all(&HgvsVariant::Rna(rna))
            .expect_err("r. input without parent must error");
        match err {
            FerroError::UnsupportedProjection { reason } => {
                assert!(reason.contains("r."), "unexpected error message: {reason}");
            }
            other => panic!("expected UnsupportedProjection, got {other:?}"),
        }
    }

    /// `HgvsVariant::Allele` whose first inner variant is c. — the
    /// fan-out path should unwrap and process the inner correctly,
    /// matching the behavior of g.-inner alleles.
    #[test]
    fn project_variant_all_accepts_allele_of_coding_input() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, HgvsVariant, LocEdit};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![HgvsVariant::Cds(cds)]));
        let results = vp
            .project_variant_all(&allele)
            .expect("Allele wrapping a c. inner must reach the fan-out path");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected the input transcript's projection from the allele's inner c., got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    #[test]
    fn project_all_no_overlap_returns_empty() {
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // g.5000 is far outside all transcripts.
        let results = vp
            .project_all("NC_000001.11:g.5000A>G")
            .expect("project_all should return Ok for no overlaps");

        assert!(
            results.is_empty(),
            "expected empty result for non-overlapping variant"
        );
    }

    // -------------------------------------------------------------------------
    // Allele compound projection tests
    // -------------------------------------------------------------------------

    /// Helper: build a cis allele `[chr1:g.1003C>A;chr1:g.1006T>A]` and
    /// project it onto NM_TEST.1.
    fn project_cis_allele() -> VariantProjection {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // Use parse → project_variant to build the allele naturally.
        // g.1003 is C (c.4) and g.1006 is T (c.7) in the test fixture
        // "ATGCGCTAA"; using the correct ref bases keeps the test
        // valid HGVS rather than just exercising control flow.
        let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("v1 parse");
        let v2 = crate::parse_hgvs("chr1:g.1006T>A").expect("v2 parse");
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v1, v2]));
        vp.project_variant(&allele, "NM_TEST.1")
            .expect("cis allele projection should succeed")
    }

    #[test]
    fn project_cis_allele_produces_cis_coding_allele() {
        let result = project_cis_allele();
        let coding = result.coding.as_ref().expect("c. allele should be present");
        // The coding variant should itself be an Allele.
        assert!(
            matches!(coding, HgvsVariant::Allele(av) if av.phase == AllelePhase::Cis),
            "expected Cis coding allele, got: {}",
            coding
        );
        // Two inner c. variants.
        if let HgvsVariant::Allele(av) = coding {
            assert_eq!(
                av.variants.len(),
                2,
                "expected 2 inner c. variants in cis allele"
            );
        }
    }

    #[test]
    fn project_cis_allele_aggregates_noncoding_axis() {
        // c↔n axis: every inner member carries an n. form, so the allele
        // projection exposes a cis n. allele mirroring the c. allele.
        let result = project_cis_allele();
        let noncoding = result
            .noncoding
            .as_ref()
            .expect("n. allele should be aggregated");
        assert!(
            matches!(noncoding, HgvsVariant::Allele(av) if av.phase == AllelePhase::Cis),
            "expected Cis noncoding allele, got: {noncoding}"
        );
        if let HgvsVariant::Allele(av) = noncoding {
            assert_eq!(
                av.variants.len(),
                2,
                "expected 2 inner n. variants in cis allele"
            );
            assert!(
                av.variants.iter().all(|v| matches!(v, HgvsVariant::Tx(_))),
                "every inner member of the n. allele should be an n. (Tx) variant"
            );
            // The fixture's inner members are g.1003→c.4 and g.1006→c.7 on
            // NM_TEST.1 (CDS at tx pos 1), so the aggregated n. coordinates are
            // n.4 and n.7 — pin them so an off-by-one in the derivation fails.
            let rendered: Vec<String> = av.variants.iter().map(|v| v.to_string()).collect();
            assert!(
                rendered.iter().any(|s| s.contains(":n.4")),
                "expected an n.4 member, got {rendered:?}"
            );
            assert!(
                rendered.iter().any(|s| s.contains(":n.7")),
                "expected an n.7 member, got {rendered:?}"
            );
        }
    }

    #[test]
    fn project_trans_allele_preserves_trans_phase() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // Ref bases match the fixture: g.1003 = C, g.1006 = T.
        let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("v1 parse");
        let v2 = crate::parse_hgvs("chr1:g.1006T>A").expect("v2 parse");
        let allele = HgvsVariant::Allele(AlleleVariant::trans(vec![v1, v2]));
        let result = vp
            .project_variant(&allele, "NM_TEST.1")
            .expect("trans allele projection should succeed");

        let coding = result.coding.as_ref().expect("c. expected");
        assert!(
            matches!(coding, HgvsVariant::Allele(av) if av.phase == AllelePhase::Trans),
            "expected Trans coding allele, got: {}",
            coding
        );
    }

    #[test]
    fn project_allele_with_frameshift_inner_sets_is_frameshift() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // g.1004del is a 1-base deletion → frameshift.
        let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("v1 parse");
        let v_fs = crate::parse_hgvs("NC_000001.11:g.1004del").expect("fs parse");
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v1, v_fs]));
        let result = vp
            .project_variant(&allele, "NM_TEST.1")
            .expect("allele with frameshift should project");

        assert!(
            result.is_frameshift,
            "allele containing a frameshift inner variant should set is_frameshift=true"
        );
    }

    #[test]
    fn project_allele_with_non_protein_inner_has_no_protein() {
        // An allele where one inner variant has no protein (intronic) →
        // the whole allele protein should be None.
        let (projector, provider) = make_intronic_test_data();
        let vp = VariantProjector::new(projector, provider);

        // g.1003 (exonic, ref C in NM_INTR.1 fixture) + g.1015 (intronic, ref N
        // since g.1015 is in the intron gap — the projector doesn't validate
        // intronic ref bases against the genomic reference).
        let v_exon = crate::parse_hgvs("NC_000001.11:g.1003C>G").expect("exon parse");
        let v_intron = crate::parse_hgvs("NC_000001.11:g.1015N>G").expect("intron parse");
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v_exon, v_intron]));
        let result = vp
            .project_variant(&allele, "NM_INTR.1")
            .expect("allele with intronic variant should project");

        assert!(
            result.protein.is_none(),
            "allele with intronic inner variant should have no protein"
        );
        assert!(result.is_intronic, "should be marked intronic");
    }

    // -------------------------------------------------------------------------
    // project_to_genomic tests (#327)
    // -------------------------------------------------------------------------

    mod project_to_genomic_tests {
        use super::*;
        use crate::hgvs::edit::{Base, InsertedSequence, NaEdit};
        use crate::hgvs::interval::{CdsInterval, GenomeInterval, TxInterval};
        use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
        use crate::hgvs::variant::{
            Accession, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, TxVariant,
        };

        /// Build an `NG_*` parent `Accession`.
        fn ng_parent(num: &str, version: u32) -> Accession {
            Accession::new("NG", num, Some(version))
        }

        /// Attach a genomic_context to an existing `CdsVariant`.
        fn attach_genomic_context_cds(mut v: CdsVariant, ctx: Accession) -> CdsVariant {
            v.accession = v.accession.with_genomic_context(ctx);
            v
        }

        /// Attach a genomic_context to an existing `TxVariant` (n.).
        fn attach_genomic_context_tx(mut v: TxVariant, ctx: Accession) -> TxVariant {
            v.accession = v.accession.with_genomic_context(ctx);
            v
        }

        /// Build a non-coding transcript test fixture for `n.` projection.
        ///
        /// `NR_TEST.1` is a 9-base non-coding RNA on chr1, plus strand,
        /// genome [1000, 1009).  Without a CDS the transcript uses pure n.
        /// coordinates (1-based: n.1..n.9). n.5 maps to genome 1004.
        fn make_noncoding_test_data() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NR_TEST.1".to_string(),
                CdotTranscript {
                    gene_name: Some("TESTGENE".to_string()),
                    contig: "chr1".to_string(),
                    strand: ProvStrand::Plus,
                    exons: vec![[1000, 1009, 0, 9]],
                    // No CDS for n. transcripts.
                    cds_start: None,
                    cds_end: None,
                    gene_id: None,
                    protein: None,
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);

            let mut provider = MockProvider::new();
            provider.add_transcript(Transcript {
                id: "NR_TEST.1".to_string(),
                gene_symbol: Some("TESTGENE".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAA".to_string()),
                cds_start: None,
                cds_end: None,
                exons: vec![Exon::with_genomic(1, 1, 9, 1000, 1008)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1008),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
                protein_id: None,
            });
            let prefix = "N".repeat(999);
            let suffix = "N".repeat(100);
            provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
            (projector, provider)
        }

        // Note: the non-coding (`else`) branch of `project_single_inner` — where
        // `noncoding` mirrors the `coding` `Tx` form — is currently unreachable
        // via full projection: the downstream g.→c. path errors with
        // `TranscriptNotOverlapping` for non-coding transcripts, pinned in
        // `tests/issue_328_projector_accepts_cnr.rs`. When that path is expanded
        // to support non-coding tx, add a `noncoding == coding` regression test
        // here. The branch stays for forward-compatibility (correct semantics
        // the moment the path is reachable).

        // -- 1: plus-strand substitution ---------------------------------------

        #[test]
        fn project_to_genomic_plus_strand_substitution() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // Build NG_TEST.1(NM_TEST.1):c.4C>A.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TEST", 1));
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("plus-strand c. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant, got: {}", out),
            };
            // Parent NG_TEST.1, plus-strand → g.1003C>A.
            assert_eq!(g.accession.to_string(), "NG_TEST.1");
            let s = g.to_string();
            assert!(
                s.starts_with("NG_TEST.1") && s.contains(":g.1003C>A"),
                "expected NG_TEST.1...:g.1003C>A, got: {}",
                s
            );
        }

        // -- 2: minus-strand substitution --------------------------------------

        #[test]
        fn project_to_genomic_minus_strand_substitution_revcomps() {
            let (projector, provider) = make_minus_strand_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.4C>T on minus-strand NM_TEST_MINUS.1 → tx_pos 4 on minus,
            // c.1 is the last base of the exon [1000, 1009) at genome 1008,
            // so c.4 = genome 1005.  Ref C on transcript = G on genome;
            // alt T on transcript = A on genome → g.1005G>A.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST_MINUS.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::T,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TESTMINUS", 1));
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("minus-strand c. → g. should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NG_TESTMINUS.1");
            let s = g.to_string();
            assert!(
                s.starts_with("NG_TESTMINUS.1") && s.contains(":g.1005G>A"),
                "minus-strand c.4C>T should revcomp to NG_TESTMINUS.1...:g.1005G>A, got: {}",
                s
            );
        }

        // -- 3: minus-strand insertion (revcomp inserted seq) ------------------

        #[test]
        fn project_to_genomic_minus_strand_insertion_revcomps_inserted_seq() {
            let (projector, provider) = make_minus_strand_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.4_5insATC on minus strand.  c.4 = g.1005, c.5 = g.1004 (minus).
            // Inserted "ATC" on the transcript → revcomp "GAT" on the genome.
            // g. start/end are min/max of the two genome coords, so 1004_1005.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST_MINUS.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::new(CdsPos::new(4), CdsPos::new(5)),
                    NaEdit::Insertion {
                        sequence: InsertedSequence::Literal("ATC".parse().unwrap()),
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TESTMINUS", 1));
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("minus-strand insertion should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            let s = g.to_string();
            // Pin the exact interval — a start/end swap regression would
            // otherwise still pass via `contains("ins")` + `contains("GAT")`.
            // c.4 = g.1005, c.5 = g.1004; canonical g. interval is low→high.
            assert!(
                s.contains(":g.1004_1005ins"),
                "expected interval :g.1004_1005ins, got: {}",
                s
            );
            assert!(
                s.contains("GAT"),
                "expected revcomp(ATC) = GAT in g., got: {}",
                s
            );
        }

        // -- 4: intronic offset, plus strand -----------------------------------

        #[test]
        fn project_to_genomic_intronic_offset_plus() {
            let (projector, provider) = make_intronic_test_data();
            let vp = VariantProjector::new(projector, provider);

            // NM_INTR.1 layout (plus strand):
            //   Exon 1: genome [1000, 1010), tx [0, 10)
            //   Intron: genome [1010, 2000)
            //   Exon 2: genome [2000, 2010), tx [10, 20)
            // c.10 is the last base of exon 1 (tx_pos = 9, 0-based, genome 1009).
            // c.10+5 → 5 bases into the intron from the 5' boundary → g.1014.
            let cds = CdsVariant {
                accession: parse_accession("NM_INTR.1"),
                gene_symbol: Some("INTRGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::with_offset(10, 5)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("INTR", 1));
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("intronic c.10+5del should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            // Expect a deletion at g.1014 on NG_INTR.1.
            assert_eq!(g.accession.to_string(), "NG_INTR.1");
            let start = g
                .loc_edit
                .location
                .start
                .inner()
                .expect("start should be concrete");
            assert_eq!(
                start.base, 1014,
                "expected g.1014 for c.10+5 on plus strand, got: {}",
                start.base
            );
        }

        // -- 5: idempotence on Genome input ------------------------------------

        #[test]
        fn project_to_genomic_idempotent_on_genome_input() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // Build a plain Genome variant.
            let g = GenomeVariant {
                accession: parse_accession("NC_000001.11"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(1003)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Genome(g.clone());
            let out = vp
                .project_to_genomic(&input)
                .expect("Genome input should pass through");
            assert_eq!(
                out.to_string(),
                input.to_string(),
                "Genome input should be idempotent under project_to_genomic"
            );
        }

        // -- 6: missing parent NG/NC -------------------------------------------

        #[test]
        fn project_to_genomic_missing_parent_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // NM_TEST.1:c.4C>A — no genomic_context attached.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("missing parent should error");
            match err {
                FerroError::UnsupportedProjection { reason } => {
                    assert!(
                        reason.contains("parent reference"),
                        "expected reason to mention 'parent reference', got: {}",
                        reason
                    );
                }
                other => panic!("expected UnsupportedProjection, got: {:?}", other),
            }
        }

        // -- 7: n. (non-coding) input ------------------------------------------

        #[test]
        fn project_to_genomic_tx_input_on_noncoding() {
            let (projector, provider) = make_noncoding_test_data();
            let vp = VariantProjector::new(projector, provider);

            // NR_TEST.1:n.5A>G on plus strand → n.5 = genome 1004.
            let tx = TxVariant {
                accession: parse_accession("NR_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    TxInterval::point(TxPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::G,
                    },
                ),
            };
            let tx = attach_genomic_context_tx(tx, ng_parent("NRTEST", 1));
            let input = HgvsVariant::Tx(tx);
            let out = vp
                .project_to_genomic(&input)
                .expect("n. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NG_NRTEST.1");
            let s = g.to_string();
            assert!(
                s.starts_with("NG_NRTEST.1") && s.contains(":g.1004A>G"),
                "expected NG_NRTEST.1...:g.1004A>G, got: {}",
                s
            );
        }

        // -- 8: r. input -------------------------------------------------------

        #[test]
        fn project_to_genomic_rna_input_with_dna_bases() {
            use crate::hgvs::interval::RnaInterval;
            use crate::hgvs::location::RnaPos;
            use crate::hgvs::variant::RnaVariant;

            // Verifies r. → g. coordinate projection on the plus strand with
            // DNA-typed bases (the case ferro emits internally after parsing).
            // RnaPos shape mirrors CdsPos so the coordinate math is identical.
            //
            // Plus-strand `Base::U` in r. inputs is NOT auto-translated to
            // `Base::T` on the g. output today — the impl forwards the edit
            // verbatim on plus strand. Translating U→T at the r. boundary is
            // a follow-up (tracked separately); for now the impl assumes the
            // caller supplies DNA-typed bases. Test below pins the working
            // path and intentionally leaves the U-input gap uncovered.
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // r.4 on NM_TEST.1 (plus strand) → g.1003 on NG_TEST.1.
            let rna = RnaVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    RnaInterval::point(RnaPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let mut rna = rna;
            rna.accession = rna.accession.with_genomic_context(ng_parent("TEST", 1));
            let input = HgvsVariant::Rna(rna);
            let out = vp
                .project_to_genomic(&input)
                .expect("r. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            // Pin both the interval and the bases: a coordinate regression
            // would otherwise pass `contains("C>A")` alone.
            let s = g.to_string();
            assert!(
                s.contains(":g.1003C>A"),
                "expected NG_TEST.1...:g.1003C>A, got: {}",
                s
            );
            assert_eq!(g.accession.to_string(), "NG_TEST.1");
        }

        // -- 9: protein input → unsupported ------------------------------------

        #[test]
        fn project_to_genomic_unsupported_for_protein() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            let input = crate::parse_hgvs("NP_TEST.1:p.Arg2Cys").expect("p. parse should succeed");
            assert!(matches!(input, HgvsVariant::Protein(_)));
            let err = vp
                .project_to_genomic(&input)
                .expect_err("p. → g. should be unsupported");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected UnsupportedProjection for p., got: {:?}",
                err
            );
        }

        // -- 10: allele input → unsupported ------------------------------------

        #[test]
        fn project_to_genomic_unsupported_for_allele() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("parse v1");
            let v2 = crate::parse_hgvs("chr1:g.1006T>A").expect("parse v2");
            let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v1, v2]));
            let err = vp
                .project_to_genomic(&allele)
                .expect_err("allele inputs should be rejected (see #328)");
            match err {
                FerroError::UnsupportedProjection { reason } => {
                    assert!(
                        reason.to_ascii_lowercase().contains("allele"),
                        "expected reason to mention 'allele', got: {}",
                        reason
                    );
                }
                other => panic!("expected UnsupportedProjection, got: {:?}", other),
            }
        }

        // -- 11: ? position sentinel → unsupported -----------------------------

        #[test]
        fn project_to_genomic_unknown_position_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // Build NG_TEST.1(NM_TEST.1):c.?C>A — start/end are CdsPos::unknown(None).
            let unknown_pos = CdsPos::unknown(None);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1").with_genomic_context(ng_parent("TEST", 1)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(unknown_pos),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("? position should be unsupported");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected UnsupportedProjection for unknown position, got: {:?}",
                err
            );
        }

        // -- 12: parsed `c.?` routes to UnsupportedProjection ------------------

        /// Pins that a *parsed* unknown-position input (`UncertainBoundary::
        /// Single(Mu::Unknown)`) routes to `UnsupportedProjection`, not
        /// `InvalidCoordinates`. Test 11 above covers the synthetic
        /// `CdsPos::unknown()` shape; this one covers the parser-produced
        /// shape, which has different AST internals and previously slipped
        /// through `.inner()` as `None` indistinguishable from a Range.
        #[test]
        fn project_to_genomic_parsed_question_position_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // `c.?_10del` parses to a Cds variant whose start endpoint is
            // `UncertainBoundary::Single(Mu::Unknown)`. We need to attach a
            // genomic_context so the variant kind is correct; reading the
            // parsed variant out and mutating its accession is sufficient.
            let parsed = crate::parse_hgvs("NM_TEST.1:c.?_10del").expect("parse should succeed");
            let cds = match parsed {
                HgvsVariant::Cds(c) => c,
                other => panic!("expected Cds variant, got: {:?}", other),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TEST", 1));
            let input = HgvsVariant::Cds(cds);

            let err = vp
                .project_to_genomic(&input)
                .expect_err("parsed c.? should be unsupported");
            match err {
                FerroError::UnsupportedProjection { reason } => {
                    assert!(
                        reason.contains("unknown position") || reason.contains("`?`"),
                        "expected reason to mention unknown position, got: {}",
                        reason
                    );
                }
                other => panic!(
                    "expected UnsupportedProjection for parsed c.?_10del, got: {:?}",
                    other
                ),
            }
        }

        // -- 13 + 14: UTR position projection (c.*N, c.-N) ---------------------

        /// Build a fixture with 5'UTR (3 bases) + CDS (9 bases) + 3'UTR (3 bases) so
        /// c.-1, c.1, c.9, c.*1 etc. all have well-defined genomic counterparts.
        ///
        /// Layout (plus strand):
        ///   Exon 1: genome [1000, 1015), tx [0, 15)  — 15 bp
        ///   cds_start = 3 (0-based)  ⇒ c.-3 = g.1000, c.-1 = g.1002
        ///   cds_end   = 12 (exclusive) ⇒ c.1 = g.1003, c.9 = g.1011, c.*1 = g.1012
        ///   sequence  = "TTTATGCGCTAAGGG" (5'UTR + CDS + 3'UTR)
        fn make_utr_test_provider_and_projector() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NM_UTR.1".to_string(),
                CdotTranscript {
                    gene_name: Some("UTRGENE".to_string()),
                    contig: "chr1".to_string(),
                    strand: ProvStrand::Plus,
                    exons: vec![[1000, 1015, 0, 15]],
                    cds_start: Some(3),
                    cds_end: Some(12),
                    gene_id: None,
                    protein: Some("NP_UTR.1".to_string()),
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);

            let mut provider = MockProvider::new();
            provider.add_transcript(Transcript {
                id: "NM_UTR.1".to_string(),
                gene_symbol: Some("UTRGENE".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("TTTATGCGCTAAGGG".to_string()),
                cds_start: Some(4), // 1-based inclusive
                cds_end: Some(12),
                exons: vec![Exon::with_genomic(1, 1, 15, 1000, 1014)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1014),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
                protein_id: Some("NP_UTR.1".to_string()),
            });
            let prefix = "N".repeat(1000);
            let suffix = "N".repeat(100);
            provider.add_genomic_sequence("chr1", format!("{}TTTATGCGCTAAGGG{}", prefix, suffix));
            (projector, provider)
        }

        #[test]
        fn project_to_genomic_3utr_star_position() {
            let (projector, provider) = make_utr_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.*1A>G on NM_UTR.1 → tx_pos = cds_end (12) → genome 1012.
            // (Without the utr3-aware fix in data::mapping::cds_to_genome, ferro
            // silently mapped c.*1 as c.1 and produced g.1003.)
            let cds = CdsVariant {
                accession: parse_accession("NM_UTR.1"),
                gene_symbol: Some("UTRGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 1,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::G,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("UTR", 1));
            let input = HgvsVariant::Cds(cds);

            let out = vp
                .project_to_genomic(&input)
                .expect("c.*1A>G should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NG_UTR.1");
            let s = g.to_string();
            assert!(
                s.contains(":g.1012A>G"),
                "expected NG_UTR.1...:g.1012A>G for c.*1A>G, got: {}",
                s
            );
        }

        #[test]
        fn project_to_genomic_5utr_negative_position() {
            let (projector, provider) = make_utr_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.-1T>G on NM_UTR.1 → tx_pos = cds_start - 1 = 2 → genome 1002
            // (seq[2] = 'T' in "TTTATGCGCTAAGGG").
            let cds = CdsVariant {
                accession: parse_accession("NM_UTR.1"),
                gene_symbol: Some("UTRGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(-1)),
                    NaEdit::Substitution {
                        reference: Base::T,
                        alternative: Base::G,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("UTR", 1));
            let input = HgvsVariant::Cds(cds);

            let out = vp
                .project_to_genomic(&input)
                .expect("c.-1T>G should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NG_UTR.1");
            let s = g.to_string();
            assert!(
                s.contains(":g.1002T>G"),
                "expected NG_UTR.1...:g.1002T>G for c.-1T>G, got: {}",
                s
            );
        }

        // -- 15: intronic offset on non-coding transcript → unsupported (#332) -

        // -- special positions (pter/qter/cen) in project_to_genomic ----------

        /// A `pter` CDS position is special (base==0, special.is_some()): it must
        /// not slip past the `is_unknown()` guard and then fail with "Invalid CDS
        /// position: 0". The extended guard must catch it as `UnsupportedProjection`.
        #[test]
        fn project_to_genomic_special_pter_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1").with_genomic_context(ng_parent("TEST", 1)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::pter()),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("pter position must not project to genomic");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected UnsupportedProjection for pter position, got: {:?}",
                err
            );
        }

        /// A `qter` endpoint must be rejected with the same error.
        #[test]
        fn project_to_genomic_special_qter_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1").with_genomic_context(ng_parent("TEST", 1)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::new(CdsPos::new(1), CdsPos::qter()),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("qter position must not project to genomic");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected UnsupportedProjection for qter position, got: {:?}",
                err
            );
        }

        /// Pins the deferred-to-#332 limitation: intronic offsets on `n.`
        /// (non-coding) inputs return `UnsupportedProjection` with a reason
        /// that references #332. If/when #332 lands and the runner is wired
        /// up, this test should start failing (XPASS) and be promoted.
        #[test]
        fn project_to_genomic_intronic_offset_on_noncoding_unsupported() {
            let (projector, provider) = make_noncoding_test_data();
            let vp = VariantProjector::new(projector, provider);

            // n.5+3del on NR_TEST.1 — non-coding transcript with intronic offset.
            let tx = TxVariant {
                accession: parse_accession("NR_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    TxInterval::point(TxPos::with_offset(5, 3)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let tx = attach_genomic_context_tx(tx, ng_parent("NRTEST", 1));
            let input = HgvsVariant::Tx(tx);

            let err = vp
                .project_to_genomic(&input)
                .expect_err("intronic n. should be unsupported pending #332");
            match err {
                FerroError::UnsupportedProjection { reason } => {
                    assert!(
                        reason.contains("#332"),
                        "expected reason to reference #332, got: {}",
                        reason
                    );
                }
                other => panic!("expected UnsupportedProjection, got: {:?}", other),
            }
        }
    }

    // -------------------------------------------------------------------------
    // Multi-build cdot build-hint propagation (closes #389 item 2)
    //
    // Both `project_to_genomic` (c./n./r. → g.) and `project_single_inner`
    // (g. → c./n./r.) MUST consult the genome build implied by the
    // input's NG/NC parent (resp. the g. variant's NC_* accession) when
    // looking up the transcript and converting coordinates. Pre-#389
    // those sites used `cdot.get_transcript` / `transcript_genome_span` /
    // `cds_to_genome` / `genome_to_cds` directly, which always returned
    // the primary build's view — silently mis-projecting GRCh37 inputs
    // against a GRCh38-primary cdot (and vice versa).
    // -------------------------------------------------------------------------
    mod multi_build_projection_tests {
        use super::*;
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::{CdsInterval, GenomeInterval};
        use crate::hgvs::location::{CdsPos, GenomePos};
        use crate::hgvs::variant::{Accession, CdsVariant, GenomeVariant, HgvsVariant, LocEdit};

        /// Minimal multi-build cdot fixture: one coding transcript whose
        /// GRCh37 view sits at genome [10000, 10010) on NC_000001.10 and
        /// whose GRCh38 view sits at genome [20000, 20010) on
        /// NC_000001.11. Different genomic coords across builds is what
        /// lets us tell which view the projector consulted.
        ///
        /// `start_codon = 0`, `stop_codon = 10` ⇒ the whole exon is CDS
        /// (no UTRs), so c.1..c.10 maps to tx_pos 0..9. The fixture is
        /// coding so `project_single_inner`'s `genome_to_cds_on_build`
        /// call lands on the CDS branch — non-coding would hit a separate
        /// "Cannot convert tx position to CDS" error unrelated to #389.
        fn multi_build_cdot_json() -> &'static str {
            r#"
            {
                "transcripts": {
                    "NM_MB_TEST.1": {
                        "gene_name": "MBTEST",
                        "genome_builds": {
                            "GRCh37": {
                                "contig": "NC_000001.10",
                                "strand": "+",
                                "exons": [[10000, 10010, 1, 0, 10, "M10"]]
                            },
                            "GRCh38": {
                                "contig": "NC_000001.11",
                                "strand": "+",
                                "exons": [[20000, 20010, 1, 0, 10, "M10"]]
                            }
                        },
                        "start_codon": 0,
                        "stop_codon": 10
                    }
                }
            }
            "#
        }

        /// Build a `VariantProjector` whose cdot mapper has the two-build
        /// fixture above. `primary` selects which build the primary
        /// transcripts table is populated from (the other build lives in
        /// `alt_build_transcripts`). The `MockProvider` carries
        /// NM_MB_TEST.1 so that downstream protein prediction inside
        /// `project_single_inner` can fetch the transcript without
        /// hitting `ReferenceNotFound`.
        fn make_multi_build_projector(primary: &str) -> VariantProjector<MockProvider> {
            use crate::reference::transcript::{
                Exon, GenomeBuild as RefGenomeBuild, ManeStatus, Strand as TxStrand,
                Transcript as RefTranscript,
            };
            use std::sync::OnceLock;

            let cdot =
                CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), primary)
                    .expect("multi-build cdot JSON should parse");
            let projector = Projector::new(cdot);
            let mut provider = MockProvider::new();
            provider.add_transcript(RefTranscript {
                id: "NM_MB_TEST.1".to_string(),
                gene_symbol: Some("MBTEST".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAATC".to_string()),
                cds_start: Some(1),
                cds_end: Some(10),
                exons: vec![Exon::new(1, 1, 10)],
                chromosome: Some("NC_000001.11".to_string()),
                genomic_start: Some(20000),
                genomic_end: Some(20009),
                genome_build: RefGenomeBuild::GRCh38,
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
            VariantProjector::new(projector, provider)
        }

        fn nc_chr1(version: u32) -> Accession {
            Accession::new("NC", "000001", Some(version))
        }

        fn ng_parent(num: &str, version: u32) -> Accession {
            Accession::new("NG", num, Some(version))
        }

        /// Smoke test the JSON parses into a usable multi-build mapper —
        /// keeps the later test failures readable if the JSON ever drifts.
        #[test]
        fn fixture_loads_both_builds() {
            let cdot =
                CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), "GRCh38")
                    .expect("multi-build JSON must parse");
            assert_eq!(cdot.primary_build(), Some("GRCh38"));
            assert_eq!(cdot.transcript_count(), 1);
            let tx38 = cdot
                .get_transcript_on_build("NM_MB_TEST.1", "GRCh38")
                .expect("GRCh38 view must be populated in primary transcripts");
            assert_eq!(tx38.contig, "NC_000001.11");
            assert_eq!(tx38.exons, vec![[20000, 20010, 0, 10]]);
            assert_eq!(tx38.cds_start, Some(0));
            assert_eq!(tx38.cds_end, Some(10));
            let tx37 = cdot
                .get_transcript_on_build("NM_MB_TEST.1", "GRCh37")
                .expect("GRCh37 view must be populated in alt_build_transcripts");
            assert_eq!(tx37.contig, "NC_000001.10");
            assert_eq!(tx37.exons, vec![[10000, 10010, 0, 10]]);
        }

        /// c.5A>T against NC_000001.10(NM_MB_TEST.1) → expect GRCh37
        /// coordinates (g.10004; cdot's internal exon[0] is HGVS-numeric,
        /// so c.1 sits at g.10000 and c.5 at g.10004). Pre-fix, this
        /// would have used the primary (GRCh38) view and produced
        /// ~g.20004, despite the parent's GRCh37 version.
        #[test]
        fn project_to_genomic_uses_grch37_view_for_nc_v10_parent() {
            let vp = make_multi_build_projector("GRCh38");
            let cds = CdsVariant {
                accession: parse_accession("NM_MB_TEST.1").with_genomic_context(nc_chr1(10)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("GRCh37 parent should project against GRCh37 exon table");
            match out {
                HgvsVariant::Genome(g) => {
                    let s = g.to_string();
                    assert!(
                        s.contains(":g.10004A>T"),
                        "expected GRCh37 coord g.10004A>T, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// c.5A>T against NC_000001.11(NM_MB_TEST.1) → expect GRCh38
        /// coordinates (g.20004). Sibling to the GRCh37 case; together
        /// they prove the parent's build version is what drives the
        /// cdot view selection.
        #[test]
        fn project_to_genomic_uses_grch38_view_for_nc_v11_parent() {
            let vp = make_multi_build_projector("GRCh38");
            let cds = CdsVariant {
                accession: parse_accession("NM_MB_TEST.1").with_genomic_context(nc_chr1(11)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("GRCh38 parent should project against GRCh38 exon table");
            match out {
                HgvsVariant::Genome(g) => {
                    let s = g.to_string();
                    assert!(
                        s.contains(":g.20004A>T"),
                        "expected GRCh38 coord g.20004A>T, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// g. → c.: NC_000001.10:g.10004A>T projected onto NM_MB_TEST.1
        /// → expect c.5A>T. Pre-fix, `project_single_inner` used the
        /// primary (GRCh38) `transcript_genome_span` of [20000, 20010),
        /// so position 10004 fell outside and the call returned
        /// `TranscriptNotOverlapping`. Post-fix the GRCh37 span
        /// [10000, 10010) is consulted and the projection succeeds.
        #[test]
        fn project_normalized_uses_grch37_view_for_nc_v10_input() {
            let vp = make_multi_build_projector("GRCh38");
            let g = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10004)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projection = vp
                .project_normalized(&HgvsVariant::Genome(g), "NM_MB_TEST.1")
                .expect("GRCh37 g. input must project via the alt-build cdot view");
            let coding = projection
                .coding
                .as_ref()
                .expect("coding (c.) form of the projection must be present");
            let s = coding.to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T for NC_000001.10 g. input, got: {}",
                s
            );
        }

        /// NG_* parents are build-agnostic per the HGVS spec — the
        /// projector falls back to the primary cdot view rather than
        /// forcing an alt-build lookup. With primary=GRCh38 the c.5
        /// position must therefore come back at the GRCh38 exon's
        /// position (g.20004) regardless of any GRCh37 alt-build data.
        #[test]
        fn project_to_genomic_with_ng_parent_uses_primary_build() {
            let vp = make_multi_build_projector("GRCh38");
            let cds = CdsVariant {
                accession: parse_accession("NM_MB_TEST.1")
                    .with_genomic_context(ng_parent("MBTEST", 1)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("NG-parented variant should project via the primary cdot view");
            match out {
                HgvsVariant::Genome(g) => {
                    let s = g.to_string();
                    assert!(
                        s.starts_with("NG_MBTEST.1"),
                        "expected output to inherit NG_ parent accession, got: {}",
                        s
                    );
                    assert!(
                        s.contains(":g.20004A>T"),
                        "NG (build-agnostic) parent must use primary (GRCh38) coords, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// `transcript_cache` and `ref_protein_cache` must partition by
        /// genome-build identity, not collapse to a single key per
        /// `transcript_id` (#389 follow-up).
        ///
        /// Pre-fix: the `coding` variant built inside `project_single_inner`
        /// dropped the genomic accession (`parse_accession(transcript_id)`
        /// carries no `genomic_context`), so `parent_key_for(&coding)`
        /// returned `None` for both `NC_*.10` and `NC_*.11` inputs. The
        /// two builds therefore shared cache entries — a GRCh38-built
        /// `Transcript` / `RefProteinBundle` could be served for a GRCh37
        /// projection and vice versa.
        ///
        /// Post-fix: the parallel `cache_variant` stamps the projected g.
        /// accession (`NC_*.10` vs `NC_*.11`) as `genomic_context`, so
        /// `parent_key_for(&cache_variant)` produces distinct keys for
        /// the two builds and the caches keep them apart. The test
        /// projects an indel (deletion) on each build so both branches
        /// of the protein-prediction dispatch
        /// (`cached_get_transcript_for_variant` + `cached_ref_translation`)
        /// execute.
        #[test]
        fn caches_partition_by_genome_build_for_nc_inputs() {
            let vp = make_multi_build_projector("GRCh38");

            // c.5del on GRCh37 (g.10004del on NC_000001.10).
            let g37 = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10004)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            vp.project_normalized(&HgvsVariant::Genome(g37), "NM_MB_TEST.1")
                .expect("GRCh37 deletion projection must succeed");

            // c.5del on GRCh38 (g.20004del on NC_000001.11).
            let g38 = GenomeVariant {
                accession: nc_chr1(11),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(20004)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            vp.project_normalized(&HgvsVariant::Genome(g38), "NM_MB_TEST.1")
                .expect("GRCh38 deletion projection must succeed");

            let tx_keys: std::collections::HashSet<_> = vp
                .transcript_cache
                .read()
                .expect("transcript cache poisoned")
                .keys()
                .cloned()
                .collect();
            assert_eq!(
                tx_keys.len(),
                2,
                "transcript_cache must keep GRCh37 and GRCh38 entries separate, \
                 got keys: {:?}",
                tx_keys
            );
            assert!(
                tx_keys.contains(&("NM_MB_TEST.1".to_string(), Some("NC_000001.10".to_string()))),
                "expected an NC_000001.10 key, got: {:?}",
                tx_keys
            );
            assert!(
                tx_keys.contains(&("NM_MB_TEST.1".to_string(), Some("NC_000001.11".to_string()))),
                "expected an NC_000001.11 key, got: {:?}",
                tx_keys
            );

            let rp_keys: std::collections::HashSet<_> = vp
                .ref_protein_cache
                .read()
                .expect("ref-protein cache poisoned")
                .keys()
                .cloned()
                .collect();
            assert_eq!(
                rp_keys.len(),
                2,
                "ref_protein_cache must keep GRCh37 and GRCh38 entries separate, \
                 got keys: {:?}",
                rp_keys
            );
        }

        /// Same as above but for NC_000001.11 (GRCh38, the primary build).
        /// Confirms the primary-build path still works alongside the new
        /// alt-build dispatch.
        #[test]
        fn project_normalized_uses_grch38_view_for_nc_v11_input() {
            let vp = make_multi_build_projector("GRCh38");
            let g = GenomeVariant {
                accession: nc_chr1(11),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(20004)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projection = vp
                .project_normalized(&HgvsVariant::Genome(g), "NM_MB_TEST.1")
                .expect("GRCh38 g. input must project via the primary cdot view");
            let coding = projection
                .coding
                .as_ref()
                .expect("coding (c.) form of the projection must be present");
            let s = coding.to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T for NC_000001.11 g. input, got: {}",
                s
            );
        }

        /// End-to-end: `project_variant_all` on a GRCh37 genome input must
        /// DISCOVER the overlapping transcript (no transcript id supplied) and
        /// project it. This is the corpus's actual path: the input is a bare
        /// `NC_*.10:g.…` HGVSg string and the projector must route overlap
        /// discovery to the GRCh37 (secondary) build via the contig, then
        /// build-route the downstream genome→cds projection to the GRCh37 exon
        /// table. Pre-fix, overlap discovery only indexed the primary (GRCh38)
        /// build, so `project_variant_all` returned an empty `Vec` ("got []").
        #[test]
        fn project_variant_all_discovers_grch37_transcript_by_overlap() {
            let vp = make_multi_build_projector("GRCh38");
            // NC_000001.10 is GRCh37; the transcript id is NOT supplied — it
            // must be discovered by overlap.
            let g = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10004)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projections = vp
                .project_variant_all(&HgvsVariant::Genome(g))
                .expect("project_variant_all must succeed for a GRCh37 g. input");

            assert!(
                !projections.is_empty(),
                "project_variant_all must discover the overlapping GRCh37 transcript \
                 (pre-fix this returned [])"
            );
            let proj = projections
                .iter()
                .find(|p| p.transcript_id == "NM_MB_TEST.1")
                .expect("NM_MB_TEST.1 must be among the discovered transcripts");
            let coding = proj
                .coding
                .as_ref()
                .expect("coding (c.) form must be present on the discovered projection");
            let s = coding.to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T from the GRCh37-discovered projection, got: {}",
                s
            );
        }

        /// Sibling to the above on the primary (GRCh38) build: overlap
        /// discovery on NC_000001.11 must keep working through the same path.
        #[test]
        fn project_variant_all_discovers_grch38_transcript_by_overlap() {
            let vp = make_multi_build_projector("GRCh38");
            let g = GenomeVariant {
                accession: nc_chr1(11),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(20004)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projections = vp
                .project_variant_all(&HgvsVariant::Genome(g))
                .expect("project_variant_all must succeed for a GRCh38 g. input");
            let proj = projections
                .iter()
                .find(|p| p.transcript_id == "NM_MB_TEST.1")
                .expect("NM_MB_TEST.1 must be discovered on the primary build too");
            let s = proj
                .coding
                .as_ref()
                .expect("coding form present")
                .to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T from the GRCh38-discovered projection, got: {}",
                s
            );
        }
    }

    #[test]
    fn project_empty_allele_unknown_transcript_returns_reference_not_found() {
        // Regression: the empty-allele fast path used to silently return Ok(...)
        // for any transcript_id, even one not present in the cdot mapper. That
        // was inconsistent with every other projection path, which surfaces
        // ReferenceNotFound for unknown accessions.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let empty_allele = HgvsVariant::Allele(AlleleVariant::cis(vec![]));
        let err = vp
            .project_variant(&empty_allele, "NM_DOES_NOT_EXIST.1")
            .expect_err("empty allele on unknown transcript should error");
        assert!(
            matches!(err, FerroError::ReferenceNotFound { ref id } if id == "NM_DOES_NOT_EXIST.1"),
            "expected ReferenceNotFound for unknown transcript, got: {:?}",
            err
        );
    }

    // ------------------------------------------------------------------
    // Direct c.→p. path for bare-NM_ inputs (#498)
    //
    // A bare transcript coding input (no genomic_context parent) carries
    // no genome alignment, so the c.→g.→CDS roundtrip cannot run. Protein
    // prediction only needs the CDS sequence + the 1-based CDS position,
    // which for an exonic c. variant the input already provides. These
    // tests pin that the direct path predicts protein and reports
    // `genomic = None` (there is no genomic form for a bare-NM_ input).
    // ------------------------------------------------------------------

    #[test]
    fn project_bare_nm_substitution_predicts_protein_without_genome() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4C>A — bare, no genomic_context parent.
        let variant = make_coding_variant("NM_TEST.1", None);
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ projection should succeed via the direct c.→p. path");
        // No alignment ingested for a bare NM_ ⇒ no genomic representation.
        assert!(
            proj.genomic.is_none(),
            "bare-NM_ projection should have genomic = None, got {:?}",
            proj.genomic
        );
        // The coding axis is the input itself.
        assert!(proj.coding.is_some(), "coding form should be present");
        // c.4C>A: codon 2 CGC(Arg) → AGC(Ser) ⇒ p.(Arg2Ser).
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2Ser)");
    }

    #[test]
    fn project_bare_nm_substitution_populates_noncoding_axis() {
        // c↔n axis: a bare-NM_ c. SNV carries an n. (transcript-relative) form
        // even with no genome alignment. The n. edit matches the c. edit; only
        // the coordinate frame reframes.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4C>A — bare, no genomic_context parent.
        let variant = make_coding_variant("NM_TEST.1", None);
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ projection should succeed");
        assert!(proj.coding.is_some(), "coding (c.) form should be present");
        let n = proj
            .noncoding
            .as_ref()
            .expect("noncoding (n.) axis should be populated for a coding transcript");
        assert!(
            matches!(n, HgvsVariant::Tx(_)),
            "noncoding form should be an n. (Tx) variant, got {n:?}"
        );
        let n_str = n.to_string();
        assert!(
            n_str.starts_with("NM_TEST.1"),
            "n. form should render on the transcript accession, got {n_str}"
        );
        // NM_TEST.1's CDS begins at transcript position 1 (sequence ATGCGCTAA,
        // no 5'UTR), so c.4 → n.4 and the edit is unchanged.
        assert!(
            n_str.contains(":n.4C>A"),
            "c.4C>A should reframe to n.4C>A (cds_start at tx pos 1), got {n_str}"
        );
    }

    #[test]
    fn project_declines_protein_when_transcript_not_version_exact() {
        // Issue #505: the protein path reads the transcript's CDS bases
        // directly. When the reference only carries a *different* version of
        // the requested transcript (a silent version-strip substitution),
        // translating those bases would emit a wrong protein attributed to the
        // requested version. The protein path must decline to predict rather
        // than lie. The non-protein axes (coding) are unaffected.
        let (projector, mut provider) = make_test_provider_and_projector();
        // Same input as `project_bare_nm_substitution_predicts_protein_without_genome`,
        // but the provider now reports NM_TEST.1 as not available at the exact
        // requested version (e.g. only a sibling version's bases are present).
        provider.mark_non_version_exact("NM_TEST.1");
        let vp = VariantProjector::new(projector, provider);
        let variant = make_coding_variant("NM_TEST.1", None);
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("projection should still succeed; only protein prediction declines");
        assert!(
            proj.protein.is_none(),
            "protein must be declined for a non-version-exact transcript, got {:?}",
            proj.protein
        );
        // The protein-only gate must not suppress the coding axis.
        assert!(
            proj.coding.is_some(),
            "coding form should still be present despite the protein gate"
        );
    }

    #[test]
    fn project_bare_nm_inframe_deletion_predicts_protein_without_genome() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4_6del — delete codon 2 (CGC = Arg) from ATG|CGC|TAA,
        // an in-frame single-residue deletion ⇒ p.(Arg2del). Bare, no parent.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(4), CdsPos::new(6)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ deletion projection should succeed via the direct path");
        assert!(
            proj.genomic.is_none(),
            "bare-NM_ projection has no genomic form"
        );
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2del)");
    }

    #[test]
    fn project_bare_nm_intronic_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4+5A>G — intronic offset on a bare NM_. Without a
        // genome alignment the intronic position cannot be placed, so
        // protein prediction is not attempted: protein = None, no panic.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::with_offset(4, 5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ intronic projection should succeed (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.protein.is_none(),
            "intronic variant should not predict a protein consequence"
        );
        assert!(proj.is_intronic, "is_intronic flag should be set");
    }

    // ------------------------------------------------------------------
    // `.genomic` is derived from inner projections, never the raw input
    // allele (#508 review). A transcript-coordinate (c./n./r.) allele has a
    // non-genomic `original`, and a bare-NM_ inner projection has no genomic
    // form at all — so the aggregated `.genomic` must come from the inner
    // projections' `.genomic` fields, dropping to `None` when any member
    // lacks one.
    // ------------------------------------------------------------------

    #[test]
    fn project_bare_nm_coding_allele_has_no_genomic() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // [c.4C>A;c.7T>A] on a bare NM_ (no genomic_context parent): each member
        // projects via the direct c.→p. path with genomic = None, so the
        // aggregated allele must report genomic = None rather than the input
        // c. allele.
        let mk = |base: i64, r: Base, a: Base| {
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(base)),
                    NaEdit::Substitution {
                        reference: r,
                        alternative: a,
                    },
                ),
            })
        };
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![
            mk(4, Base::C, Base::A),
            mk(7, Base::T, Base::A),
        ]));
        let proj = vp
            .project_variant(&allele, "NM_TEST.1")
            .expect("bare-NM_ coding allele should project via the direct path");
        assert!(
            proj.genomic.is_none(),
            "bare-NM_ coding allele must report genomic = None, got {:?}",
            proj.genomic
        );
        assert!(proj.coding.is_some(), "coding allele should be present");
    }

    #[test]
    fn project_empty_allele_reports_no_genomic() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // An empty allele has no inner variant to derive a genomic form from,
        // so `.genomic` must be None (consistent with coding/protein, which are
        // already None for the empty allele) rather than the input allele.
        let empty = HgvsVariant::Allele(AlleleVariant::cis(vec![]));
        let proj = vp
            .project_variant(&empty, "NM_TEST.1")
            .expect("empty allele on a known transcript should project");
        assert!(
            proj.genomic.is_none(),
            "empty allele must report genomic = None, got {:?}",
            proj.genomic
        );
        assert!(proj.coding.is_none(), "empty allele has no coding form");
        assert!(proj.protein.is_none(), "empty allele has no protein form");
    }

    #[test]
    fn project_bare_nm_unknown_position_is_unsupported_not_invalid() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::uncertainty::Mu;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A bare coding input whose position is the `?` sentinel
        // (`Single(Mu::Unknown)`) must classify as UnsupportedProjection — the
        // same contract `project_to_genomic` preserves — not InvalidCoordinates,
        // which is reserved for a genuinely missing/range coordinate (#508).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::with_uncertainty(Mu::Unknown, Mu::Unknown),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        // Use `project_normalized` to exercise the boundary classification in
        // `project_coding_direct` directly, without the normalizer rewriting the
        // `?` position first.
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("bare c.? must not project");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "bare c.? should be UnsupportedProjection, got {err:?}"
        );
    }

    #[test]
    fn project_bare_nm_concrete_unknown_sentinel_is_unsupported_not_utr() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A programmatically-built `CdsPos::unknown(None)` carries the `base == 0`
        // sentinel inside a `Mu::Certain` (not `Mu::Unknown`), so it slips past
        // `resolve_uncertain_boundary`. Without an explicit `is_unknown()` guard
        // the direct path would treat `base <= 0` as 5'UTR and return `Ok`,
        // diverging from `project_to_genomic`, which rejects the sentinel as
        // `UnsupportedProjection` (#508 review).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::unknown(None)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("concrete `?` sentinel must not project as UTR");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "concrete `?` sentinel should be UnsupportedProjection, got {err:?}"
        );
    }

    // -------------------------------------------------------------------------
    // Issue #504: CDS-boundary edits that reach the initiation codon
    // -------------------------------------------------------------------------

    /// Bare-`NM_` coding transcript with a 5 nt 5'UTR ("GCACC") preceding the
    /// CDS "ATGTACTAA" (Met-Tyr-Stop). The 5'UTR lets `c.-1` exist so an edit
    /// can span the 5'UTR→CDS boundary into the translation initiation codon.
    ///
    /// Full transcript sequence: "GCACCATGTACTAA" (14 nt).
    ///   tx 1-5  = 5'UTR  G C A C C   → c.-5 c.-4 c.-3 c.-2 c.-1
    ///   tx 6-8  = ATG    (initiation codon) → c.1 c.2 c.3
    ///   tx 9-14 = TACTAA (Tyr-Stop)
    /// So `c.-1`=C(tx5), `c.1`=A(tx6), `c.2`=T(tx7), `c.3`=G(tx8).
    fn make_utr5_provider_and_projector() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_UTR5.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // [genome_start(0-based), genome_end(0-based excl), tx_start(1-based), tx_end]
                exons: vec![[1000, 1014, 0, 14]],
                cds_start: Some(5), // 0-based: CDS starts at tx_pos 5 (5 nt 5'UTR)
                cds_end: Some(14),  // 0-based exclusive
                gene_id: None,
                protein: Some("NP_UTR5.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_UTR5.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("GCACCATGTACTAA".to_string()),
            cds_start: Some(6), // 1-based inclusive: first CDS base (A of ATG)
            cds_end: Some(14),
            exons: vec![Exon::new(1, 1, 14)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1013),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: Some("NP_UTR5.1".to_string()),
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "GCACCATGTACTAA", suffix));
        (projector, provider)
    }

    /// A boundary-spanning duplication whose span reaches into the initiation
    /// codon (CDS 1–3) has an unpredictable protein consequence and must be
    /// reported as `p.(Met1?)` — even though its 5' end is in the 5'UTR
    /// (`c.-1`, CDS base ≤ 0). Before #504 the coarse `is_utr` gate dropped it
    /// to "no protein predicted" (`None`). Exercised via `project_normalized`
    /// so the canonical boundary-spanning form is fed straight to the gate.
    #[test]
    fn boundary_dup_reaching_initiation_codon_is_met1_unknown() {
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.-1_2dup duplicates c.-1,c.1,c.2 — the copy is inserted between c.2
        // and c.3, splitting the ATG initiation codon.
        let v = crate::parse_hgvs("NM_UTR5.1:c.-1_2dup").expect("parse c.-1_2dup");
        let result = vp
            .project_normalized(&v, "NM_UTR5.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .map(|p| p.to_string())
            .expect("a boundary edit reaching the initiation codon must predict a protein");
        assert_eq!(p, "NP_UTR5.1:p.(Met1?)");
    }

    /// A 5'UTR edit that does **not** reach the initiation codon (both
    /// endpoints have CDS base ≤ 0) has no predictable CDS consequence, so the
    /// refined gate must still report `None` — not a spurious `p.(Met1?)`.
    /// Guards against #504 widening the gate too far.
    #[test]
    fn pure_5utr_edit_not_reaching_initiation_codon_predicts_no_protein() {
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.-3_-2del deletes two 5'UTR bases; the CDS is untouched.
        let v = crate::parse_hgvs("NM_UTR5.1:c.-3_-2del").expect("parse c.-3_-2del");
        let result = vp
            .project_normalized(&v, "NM_UTR5.1")
            .expect("projection should succeed");
        assert!(
            result.protein.is_none(),
            "pure-5'UTR edit must not predict a protein, got {:?}",
            result.protein.map(|p| p.to_string())
        );
    }

    /// A substitution in the initiation codon proper (CDS base 1–3, not UTR)
    /// also funnels through the gate's initiation-codon short-circuit, yielding
    /// `p.(Met1?)`. Confirms #504's refactor did not regress the in-CDS
    /// start-codon path (#512).
    #[test]
    fn in_cds_initiation_codon_substitution_is_met1_unknown() {
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.1A>G changes the first base of the ATG start codon.
        let v = crate::parse_hgvs("NM_UTR5.1:c.1A>G").expect("parse c.1A>G");
        let result = vp
            .project_normalized(&v, "NM_UTR5.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .map(|p| p.to_string())
            .expect("a start-codon substitution must predict a protein");
        assert_eq!(p, "NP_UTR5.1:p.(Met1?)");
    }

    // -------------------------------------------------------------------------
    // Direct c.→p. path: special positions (pter/qter/cen) must be rejected
    // as UnsupportedProjection — not silently mis-classified as 5'UTR because
    // `base == 0` satisfies `cds_start.base <= 0`.
    // -------------------------------------------------------------------------

    #[test]
    fn project_bare_nm_special_pter_is_unsupported_not_5utr() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A `pter` CDS position has `base == 0` and `special == Some(Pter)`.
        // Without the `is_special()` guard, `cds_start.base <= 0` evaluates to
        // `true` and the variant is silently classified as 5'UTR, producing a
        // wrong protein consequence. The guard must reject it as
        // `UnsupportedProjection` (#488).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::pter()),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("pter sentinel must not project as 5'UTR");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "pter sentinel should be UnsupportedProjection, got {err:?}"
        );
    }

    #[test]
    fn project_bare_nm_special_qter_is_unsupported() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // `qter` endpoint: base == 0, special == Some(Qter). Same bug path as
        // pter — must be rejected rather than silently misclassified.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::qter()),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("qter sentinel must not project as UTR");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "qter sentinel should be UnsupportedProjection, got {err:?}"
        );
    }
}
