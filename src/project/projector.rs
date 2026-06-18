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
    affects_initiation_codon, build_initiator_unknown, build_whole_protein_unknown,
    cds_has_valid_start, predict_indel_protein, predict_substitution_protein, read_cds_start_codon,
    whole_exon_deletion_span, RefProteinBundle,
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
    /// Explicit genome-build override for build-agnostic inputs (#715).
    ///
    /// A bare `NG_`/`LRG_` input carries no build, so by default it resolves to
    /// the GRCh38-preferred placement. When set (via [`Self::with_assembly`]),
    /// this build (`"GRCh37"` / `"GRCh38"`) *fills in* for such inputs — it does
    /// **not** override a build the input accession already encodes
    /// (`NC_*.10`/`.11`, `GRCh37(...)`); the accession's build stays
    /// authoritative (see [`Self::build_hint_for_variant`]).
    assembly_override: Option<&'static str>,
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
            assembly_override: None,
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

    /// Infer the genome build the *input variant itself* encodes.
    ///
    /// Examines the variant's own accession first: a g. variant on
    /// `NC_*.10` carries GRCh37, on `NC_*.11` carries GRCh38, and an
    /// assembly-style `GRCh37(chr1):g.…` carries GRCh37 directly. If the
    /// variant has no inherent build (e.g. a c. variant on a bare NM)
    /// but does have a `genomic_context` (NG/NC parent), the parent is
    /// consulted. Anything else returns `None` (e.g. a bare `NG_`/`LRG_`,
    /// which carries no build), so the caller falls back to the
    /// build-agnostic primary-cdot path (issue #389).
    fn infer_input_build(variant: &HgvsVariant) -> Option<&'static str> {
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

    /// The genome build to project `variant` under, honoring the
    /// [`Self::with_assembly`] override (#715).
    ///
    /// The build the input *encodes* is authoritative; the `--assembly` override
    /// only *fills in* when the input is build-agnostic (a bare `NG_`/`LRG_`).
    /// A contradictory override is warned and ignored at the public entry points
    /// (see [`Self::warn_assembly_conflict`]) — it never reinterprets an explicit
    /// `NC_*.10`/`.11` accession, which would emit wrong coordinates under a
    /// right-looking accession (the #655/#702 failure class).
    fn build_hint_for_variant(&self, variant: &HgvsVariant) -> Option<&'static str> {
        Self::infer_input_build(variant).or(self.assembly_override)
    }

    /// Emit a warning when an explicit `--assembly` override contradicts the
    /// build the input accession already encodes (#715). The input wins; this
    /// only flags the ignored override. Called once per public projection entry.
    fn warn_assembly_conflict(&self, variant: &HgvsVariant) {
        if let (Some(override_build), Some(input_build)) =
            (self.assembly_override, Self::infer_input_build(variant))
        {
            if override_build != input_build {
                log::warn!(
                    "ignoring --assembly {override_build}: input accession already \
                     specifies {input_build}"
                );
            }
        }
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
                let build_hint = self.build_hint_for_variant(effective);
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
                let build_hint = self.build_hint_for_variant(effective);
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
                let build_hint = self.build_hint_for_variant(effective);
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

    /// Set an explicit genome-build override for build-agnostic inputs (#715).
    ///
    /// `build` must be a canonical ferro build string (`"GRCh37"` / `"GRCh38"`);
    /// validate/normalize user-supplied names at the boundary with
    /// [`crate::liftover::aliases::normalize_assembly_name`]. The override only
    /// *fills in* a build for inputs that don't encode one (a bare `NG_`/`LRG_`);
    /// an input whose accession already specifies a build keeps it, and a
    /// contradictory override is warned and ignored (see
    /// [`Self::build_hint_for_variant`]). `None` clears the override.
    pub fn with_assembly(mut self, build: Option<&'static str>) -> Self {
        self.assembly_override = build;
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
        self.warn_assembly_conflict(variant);
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
        self.warn_assembly_conflict(variant);
        self.project_variant_inner(variant, transcript_id)
    }

    /// Parse, normalize, and project an HGVS string onto the *curated* set of
    /// overlapping transcripts (#656), returning results in clinical priority
    /// order (MANE Select first, then Plus Clinical, then canonical, then
    /// longest CDS).
    ///
    /// The result is the curated enumerated set rather than every
    /// cdot-overlapping record: superseded transcript versions are collapsed
    /// (only the highest version per base accession is kept) and predicted
    /// `XM_`/`XR_` models are dropped when a curated `NM_`/`NR_` transcript
    /// covers the same locus; predicted models are kept only when they are the
    /// sole coverage. See [`select_enumerated_transcript_ids`] for the policy.
    ///
    /// Returns an empty `Vec` when the variant overlaps no known transcripts.
    /// Individual transcript errors are logged at trace level and silently
    /// skipped so that a single bad transcript does not abort the whole call.
    pub fn project_all(&self, hgvs_string: &str) -> Result<Vec<VariantProjection>, FerroError> {
        let variant = crate::parse_hgvs(hgvs_string)?;
        self.project_variant_all(&variant)
    }

    /// Normalize and project an already-parsed g. variant onto the *curated*
    /// set of overlapping transcripts (#656).
    ///
    /// See [`project_all`] for the curated enumeration policy, ordering, and
    /// error-handling semantics.
    pub fn project_variant_all(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VariantProjection>, FerroError> {
        self.warn_assembly_conflict(variant);
        // 0. An NG_/LRG_ genomic *input* carries parent-relative coordinates that
        //    cdot cannot map; rewrite it into the chromosome (NC_) frame so the
        //    fan-out below can enumerate and project onto the overlapping
        //    transcripts (#480 inverse). Capture the build the *original* input
        //    carries before `deanchor` shadows `variant` — the de-anchored form is
        //    on an `NC_` accession whose build would be circular to read back
        //    (#653/#713); the parent placement must be selected by the input's
        //    build, matching what `deanchor` itself used.
        let input_build = self.build_hint_for_variant(variant);
        let (variant, parent) = self.deanchor_genomic_parent_input(variant);
        // 1. Normalize once in the genome frame, then fan out across the
        //    overlapping transcripts.
        let normalized = self.normalizer.normalize(&variant)?;
        let results = self.project_normalized_all_inner(&normalized)?;
        // 2. The plain (non-parent) path returns the genome-frame result directly.
        let Some(parent) = parent else {
            return Ok(results);
        };
        // 3. Parent path: the fan-out 3'-shifts the variant in the *genome*
        //    frame, so a repeat-region edit can land at a transcript coordinate
        //    two bases off from the HGVS-correct, transcript-frame-normalized
        //    position. Re-project each transcript from its own coding form (which
        //    normalizes in the transcript frame), then re-frame under the parent.
        let placement = self
            .provider
            .genomic_placement_on_build(&parent, input_build);
        let mut framed = Vec::with_capacity(results.len());
        for r in results {
            let r = match r.coding.as_ref() {
                // Re-project from the coding form to pick up transcript-frame
                // normalization. Tolerate a failure the same way the fan-out in
                // `project_normalized_all` tolerates a failing transcript — but
                // log it rather than swallow it silently. The genome-frame `r` we
                // fall back to is itself a valid projection (only a repeat-region
                // edit can be off — by up to the repeat-unit length — from the
                // transcript-frame-optimal position); emitting it is preferable to
                // dropping the transcript outright.
                Some(coding) => match self.project_variant(coding, &r.transcript_id) {
                    Ok(reframed) => reframed,
                    Err(e) => {
                        log::trace!(
                            "project_variant_all: transcript-frame re-projection failed for {}: \
                             {}; keeping genome-frame projection",
                            r.transcript_id,
                            e
                        );
                        r
                    }
                },
                None => r,
            };
            framed.push(Self::frame_projection_owned(r, &parent, placement.as_ref()));
        }
        Ok(framed)
    }

    /// Re-frame a [`VariantProjection`] produced from a de-anchored `NC_` variant
    /// back under the original `NG_`/`LRG_` parent (#480): the transcript-relative
    /// descriptions (`coding`/`noncoding`/`protein`) gain the parent as their
    /// `genomic_context` (rendering `NG_(NM_)` / `NG_(NP_)`), and the genomic
    /// description is re-anchored into the parent's own coordinate frame.
    ///
    /// The synthesized gene-symbol selector is dropped: per #121 ferro does not
    /// emit a selector that was not in the input, and the mutalyzer `NG_(NM_)`
    /// form carries none.
    fn frame_projection_owned(
        mut result: VariantProjection,
        parent: &Accession,
        placement: Option<&crate::reference::GenomicPlacement>,
    ) -> VariantProjection {
        for field in [
            &mut result.coding,
            &mut result.noncoding,
            &mut result.protein,
        ] {
            if let Some(v) = field.as_ref() {
                *field = Some(Self::relabel_under_parent(v, parent));
            }
        }
        if let (Some(HgvsVariant::Genome(gv)), Some(placement)) =
            (result.genomic.as_ref(), placement)
        {
            match Self::reanchor_genome_to_parent(gv.clone(), placement) {
                Ok(mut reanchored) => {
                    reanchored.accession = parent.clone();
                    reanchored.gene_symbol = None;
                    result.genomic = Some(HgvsVariant::Genome(reanchored));
                }
                Err(e) => {
                    // The genomic axis cannot be re-anchored into the parent
                    // frame (endpoint outside the placed span or an uncertain
                    // position). Drop the unframable genomic axis rather than
                    // stamp the parent accession on a chromosome coordinate,
                    // which would be invalid HGVS (#655). The transcript-relative
                    // coding/noncoding/protein forms re-labeled above remain
                    // valid and are preferable to dropping the whole projection
                    // (mirrors the tolerate-and-log fallback in
                    // `project_variant_all`).
                    log::trace!(
                        "frame_projection_owned: cannot re-anchor the genomic axis into the {} \
                         parent frame: {}; reporting genomic = None",
                        parent,
                        e
                    );
                    result.genomic = None;
                }
            }
        }
        result
    }

    /// Return a clone of `variant` whose top-level accession carries `parent` as
    /// its `genomic_context` and whose synthesized gene-symbol selector is
    /// cleared. Used to frame a transcript-relative description under its
    /// `NG_`/`LRG_` parent (#480).
    fn relabel_under_parent(variant: &HgvsVariant, parent: &Accession) -> HgvsVariant {
        let mut v = variant.clone();
        match &mut v {
            HgvsVariant::Cds(c) => {
                c.accession = c.accession.clone().with_genomic_context(parent.clone());
                c.gene_symbol = None;
            }
            HgvsVariant::Tx(t) => {
                t.accession = t.accession.clone().with_genomic_context(parent.clone());
                t.gene_symbol = None;
            }
            HgvsVariant::Rna(r) => {
                r.accession = r.accession.clone().with_genomic_context(parent.clone());
                r.gene_symbol = None;
            }
            HgvsVariant::Protein(p) => {
                p.accession = p.accession.clone().with_genomic_context(parent.clone());
                p.gene_symbol = None;
            }
            _ => {}
        }
        v
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
    /// The output uses the parent accession stored in the input's
    /// `Accession.genomic_context`:
    /// - an `NC_` chromosome parent is already in the genome frame and passes
    ///   through unchanged;
    /// - an `NG_` RefSeqGene or `LRG_` parent whose chromosomal placement is
    ///   known (see [`ReferenceProvider::genomic_placement`]) is re-anchored
    ///   into the parent's own frame (#480);
    /// - an `NG_`/`LRG_` parent with **no** known placement, or whose endpoint
    ///   cannot be re-anchored, is declined (see `# Errors`) rather than
    ///   emitting chromosome (`NC_`) coordinates under the parent accession,
    ///   which would be invalid HGVS (#655).
    ///
    /// Idempotent on `Genome` input.
    ///
    /// This is the genomic-axis *output*. The internal pivot used to derive
    /// coding/protein keeps the chromosome frame — see
    /// [`Self::project_to_genomic_nc`].
    ///
    /// # Errors
    ///
    /// Returns [`FerroError::UnsupportedProjection`] when an `NG_`/`LRG_` parent
    /// has no known chromosomal placement or an endpoint falls outside the placed
    /// span, and [`FerroError::InvalidCoordinates`] when an endpoint has no single
    /// resolved position (uncertain/compound boundary) (#655). See
    /// [`Self::project_to_genomic_nc`] for the rest of the error contract
    /// (`UnsupportedProjection` / `InvalidCoordinates` / `ReferenceNotFound`).
    pub fn project_to_genomic(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        self.warn_assembly_conflict(variant);
        // A `Genome` input is already in its parent frame: an `NC_` chromosome
        // variant or an `NG_`/`LRG_` genomic variant in its own coordinates.
        // `reanchor_genome_output` assumes an `NC_`-frame pivot produced by
        // `project_to_genomic_nc`, so feeding an already-parent-framed `NG_`/`LRG_`
        // variant through it would either double-transform (mistaking the parent
        // coordinate for an `NC_` coordinate) or decline outright. Pass `Genome`
        // inputs through unchanged to preserve idempotence (#702).
        if matches!(variant, HgvsVariant::Genome(_)) {
            return Ok(variant.clone());
        }

        // #537: a `c.pter`/`c.qter` telomere-flank input resolves to the parent
        // reference's own genomic terminus, not a transcript-mapped coordinate.
        // Handle it before the cdot pivot below, which cannot (and per #534
        // must not) number a transcript-flank marker.
        if let Some(result) = self.project_cds_terminus_to_parent(variant) {
            return result;
        }

        match self.project_to_genomic_nc(variant)? {
            HgvsVariant::Genome(gv) => {
                let reanchored =
                    self.reanchor_genome_output(gv, self.build_hint_for_variant(variant))?;
                // #737: re-normalize the genomic output in its own frame — the
                // transcript-space normalization is genome-5' on a minus-strand
                // transcript, so a shiftable variant must be re-shuffled to its
                // genomic-3' position here. Falls back to the un-normalized form
                // if normalization fails (e.g. no reference bases); the dropped
                // error is `trace!`-logged so a regression is observable rather
                // than silently masked (the fallback swallows *any* normalize
                // error, not only missing bases).
                Ok(self.normalizer.normalize(&reanchored).unwrap_or_else(|e| {
                    log::trace!(
                        "genomic-frame renormalization failed for {reanchored}; \
                             emitting un-normalized re-anchored form: {e}"
                    );
                    reanchored
                }))
            }
            other => Ok(other),
        }
    }

    /// Project a `c.pter`/`c.qter` (telomere-flank marker) coding input onto its
    /// parent reference's own genomic frame (#537).
    ///
    /// `pter` denotes the 5'-most position of the parent reference (`g.1`) and
    /// `qter` the 3'-most (`g.<length>`); a `pter_qter` range spans the whole
    /// reference. Unlike a numeric `c.` coordinate these markers do not number a
    /// transcript position — PR #534 correctly refuses to do so on the `c.`
    /// axis — so they map straight to the parent reference's termini with no
    /// cdot transcript mapping. The emitted `g.1` / `g.<length>` is normalized
    /// downstream (e.g. the 3' rule rolls `g.1del` through a leading homopolymer
    /// run, matching mutalyzer's `g.3del`).
    ///
    /// Returns `None` when `variant` is not a handled terminus case — not a
    /// `Cds` input, an endpoint that is not a `pter`/`qter` marker (numeric,
    /// `?`, `cen`, or a mixed marker/numeric range), or a reversed `qter_pter` —
    /// so the caller falls through to the normal transcript-mapped projection
    /// (which declines or numbers as appropriate). Returns `Some(Err(..))` when
    /// the parent reference is unresolvable or its length is unknown to the
    /// provider (e.g. the pinned parent version is absent — a reference-coverage
    /// gap, #645/#672 — not a projection bug).
    fn project_cds_terminus_to_parent(
        &self,
        variant: &HgvsVariant,
    ) -> Option<Result<HgvsVariant, FerroError>> {
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::location::{GenomePos, SpecialPosition};
        use crate::hgvs::variant::{GenomeVariant, LocEdit};

        let HgvsVariant::Cds(v) = variant else {
            return None;
        };
        // Both endpoints must be telomere markers. A numeric or `?` endpoint, or
        // a mixed marker/numeric range (e.g. `c.pter_-51`), needs the transcript
        // mapping and is left to the normal path.
        let start = resolve_uncertain_boundary(&v.loc_edit.location.start, "c.", "start").ok()?;
        let end = resolve_uncertain_boundary(&v.loc_edit.location.end, "c.", "end").ok()?;
        let (Some(start_marker), Some(end_marker)) = (start.special, end.special) else {
            return None;
        };

        // Resolve the parent genomic reference: an explicit `genomic_context`
        // parent, or the structural `LRG_<n>` parent of a bare LRG transcript
        // (#480). No parent → leave it to the normal path, which raises the
        // canonical "no parent reference" error.
        let parent = match v.accession.genomic_context.as_deref().cloned() {
            Some(p) => p,
            None => Self::lrg_genomic_parent(&v.accession)?,
        };
        let length = match self.provider.get_sequence_length(&parent.full()) {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };
        // A zero-length parent is unreachable for any real NG/LRG/NC contig, but
        // guard the endpoint math against it (#526): `qter` would build an
        // invalid 1-based `g.0` and `pter_qter` a reversed `g.1_0`. Decline so
        // the normal path raises the canonical refusal rather than emitting a
        // malformed interval.
        if length == 0 {
            return None;
        }

        // pter → 5'-most (g.1); qter → 3'-most (g.<length>); pter_qter → the
        // whole reference. `cen` and a reversed `qter_pter` are out of scope and
        // left to the normal path (which declines).
        let (g_start, g_end) = match (start_marker, end_marker) {
            (SpecialPosition::Pter, SpecialPosition::Pter) => (1, 1),
            (SpecialPosition::Qter, SpecialPosition::Qter) => (length, length),
            (SpecialPosition::Pter, SpecialPosition::Qter) => (1, length),
            _ => return None,
        };

        let interval = GenomeInterval::new(GenomePos::new(g_start), GenomePos::new(g_end));
        let gv = GenomeVariant {
            accession: parent,
            gene_symbol: None,
            loc_edit: LocEdit::with_uncertainty(interval, v.loc_edit.edit.clone()),
        };
        Some(Ok(HgvsVariant::Genome(gv)))
    }

    /// Re-anchor a chromosome-frame genomic projection into its parent's own
    /// frame for *output* (#480), or decline.
    ///
    /// `gv` is a genome variant stamped with the parent accession from the
    /// projected transcript's `genomic_context` (as produced by
    /// [`Self::project_to_genomic_nc`]). An `NC_` chromosome parent is already in
    /// the genome frame and is returned unchanged; an `NG_`/`LRG_` parent is
    /// re-anchored into its own frame via
    /// [`ReferenceProvider::genomic_placement_on_build`].
    ///
    /// `build` is the genome build the *original input* carries (from
    /// [`Self::build_hint_for_variant`]); it selects the parent's build-appropriate
    /// placement (#653/#713). `None` (a bare `NG_` with no build signal) keeps the
    /// GRCh38-preferred behavior.
    ///
    /// Shared by the public [`Self::project_to_genomic`] (which propagates the
    /// decline as an error) and the multi-axis single-transcript path in
    /// [`Self::project_single_inner`] (which degrades the genomic axis to `None`
    /// on decline, keeping the still-valid coding/protein axes — #702).
    ///
    /// # Errors
    ///
    /// Returns [`FerroError::UnsupportedProjection`] for an `NG_`/`LRG_` parent
    /// with no known placement or an endpoint outside the placed span, and
    /// [`FerroError::InvalidCoordinates`] for an uncertain/compound endpoint
    /// (#655) — rather than emit a chromosome coordinate under the parent
    /// accession (invalid HGVS).
    fn reanchor_genome_output(
        &self,
        gv: crate::hgvs::variant::GenomeVariant,
        build: Option<&str>,
    ) -> Result<HgvsVariant, FerroError> {
        let reanchored = match self
            .provider
            .genomic_placement_on_build(&gv.accession, build)
        {
            Some(placement) => Self::reanchor_genome_to_parent(gv, &placement)?,
            // No placement: an `NC_` chromosome parent is already in the genome
            // frame and passes through unchanged, but an `NG_`/`LRG_` parent
            // cannot be re-anchored (cdot carries only the transcript's `NC_`
            // alignment). Emitting the `NC_` coordinate under the `NG_`/`LRG_`
            // accession would be invalid HGVS, so decline instead (#480/#655).
            None if &*gv.accession.prefix == "NG" || gv.accession.is_lrg() => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "cannot project to the {} parent frame: no chromosomal \
                         placement is known for this NG_/LRG_ reference, so cdot's \
                         chromosome (NC_) coordinates cannot be re-anchored into the \
                         parent's own frame (#480/#655)",
                        gv.accession.transcript_accession(),
                    ),
                });
            }
            None => gv,
        };
        Ok(HgvsVariant::Genome(reanchored))
    }

    /// Project a transcript-coordinate variant onto its parent genomic reference
    /// in the **chromosome (`NC_`) frame** — the pivot for downstream
    /// genome→CDS/protein derivation, which cdot computes in chromosome
    /// coordinates. The public [`Self::project_to_genomic`] wraps this and
    /// re-anchors the *output* into an `NG_`/`LRG_` parent's own frame (#480).
    ///
    /// The output stamps the parent accession from `Accession.genomic_context`.
    ///
    /// # Errors
    ///
    /// Returns [`FerroError::UnsupportedProjection`] for:
    /// - `p.` / `m.` / `o.` / RNA-fusion / `Allele` / `NullAllele` / `UnknownAllele`
    ///   inputs (see #328 for allele support),
    /// - transcript-coordinate inputs whose `Accession.genomic_context` is
    ///   absent (no parent NG/NC reference) — except a bare LRG transcript,
    ///   whose `LRG_<n>` parent is derived structurally (#480),
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
    /// `r.` inputs carrying `Base::U` are translated to DNA on the way to the
    /// g. output (`U`→`T` on the plus strand, `U`→`A` via complement on the
    /// minus strand), so the emitted g. variant is always valid DNA (#395 item 4).
    fn project_to_genomic_nc(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
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
        // Prefer the genome build encoded by the parent's chromosomal placement.
        // An `NG_`/`LRG_` parent accession carries no build, so without this the
        // cdot pivot below would fall back to cdot's primary build and could
        // compute coordinates on a different `NC_` build than the one the
        // placement — and the downstream re-anchor / de-anchor steps that stamp
        // `placement.nc` — assume (#646). Fall back to the parent's own build tag
        // (e.g. an explicit `NC_*` parent) when there is no placement.
        // Prefer the build the *input variant* carries (an explicit `NC_*.10`/`.11`
        // accession or `genomic_context`, or a `GRCh37(...)` assembly tag) — the
        // authoritative signal of intended build (#653/#713). Only when the input
        // is build-agnostic (e.g. a bare `NG_` parent) fall back to the build
        // implied by the GRCh38-preferred placement, then to the parent's own
        // version. Deriving the build from the variant first (rather than from the
        // chosen placement) lets the placement *selection* below honor it, instead
        // of being circularly fixed to GRCh38.
        let build_hint = self
            .build_hint_for_variant(variant)
            .or_else(|| {
                self.provider.genomic_placement(&parent).and_then(|p| {
                    crate::liftover::aliases::infer_genome_build_from_accession(&p.nc)
                })
            })
            .or_else(|| crate::liftover::aliases::infer_genome_build_from_accession(&parent));
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
                // Sequence-aware correction (#644), mirroring the g.→c. path.
                // Only simple exonic positions are corrected; intronic / special
                // positions are derived by boundary arithmetic the realignment
                // doesn't model, so leave them to the naive result.
                if p.offset.is_some() || p.special.is_some() || p.utr3 {
                    return Ok(result.variant.base);
                }
                let tx_pos = match cdot_tx.cds_to_tx(p.base) {
                    Some(t) => t,
                    None => return Ok(result.variant.base),
                };
                self.correct_genome_for_exon_indel(
                    cdot_tx,
                    &transcript_id,
                    tx_pos,
                    result.variant.base,
                )
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
                let naive = cdot_tx.tx_to_genome(tx_pos_0based).ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: format!(
                            "tx position {} not in any exon of {}",
                            p.base, transcript_id
                        ),
                    }
                })?;
                // Sequence-aware correction (#644): the inbound g.→n. path runs
                // `correct_cds_for_exon_indel` regardless of coding status, so a
                // non-coding `n.`/`r.` coordinate is corrected on the way in.
                // Apply the matching outbound correction here so it round-trips.
                self.correct_genome_for_exon_indel(cdot_tx, &transcript_id, tx_pos_0based, naive)
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

        // This NC-frame result is the **pivot** used downstream (genome→CDS via
        // cdot needs chromosome coordinates). The public `project_to_genomic`
        // re-anchors it into an NG_/LRG_ parent's own frame for the genomic
        // *output*; the pivot must stay NC (#480).
        let g_variant = GenomeVariant {
            accession: parent,
            gene_symbol,
            loc_edit: LocEdit::with_uncertainty(g_interval, g_edit_mu),
        };
        Ok(HgvsVariant::Genome(g_variant))
    }

    /// Re-anchor an NC-frame genome variant into a genomic parent's own frame
    /// using its [`GenomicPlacement`] (#480): apply the affine NC→parent
    /// transform to the interval, and reverse-complement the edit when the
    /// parent runs antiparallel to the chromosome. If an endpoint falls outside
    /// the placed span the variant is returned unchanged (chromosome frame).
    fn reanchor_genome_to_parent(
        gv: crate::hgvs::variant::GenomeVariant,
        placement: &crate::reference::GenomicPlacement,
    ) -> Result<crate::hgvs::variant::GenomeVariant, FerroError> {
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::variant::GenomeVariant;
        let bounds = match (
            gv.loc_edit.location.start.inner(),
            gv.loc_edit.location.end.inner(),
        ) {
            (Some(s), Some(e)) => Some((s.base, e.base)),
            _ => None,
        };
        let Some((nc_lo, nc_hi)) = bounds else {
            // Uncertain/compound endpoint: no single position to re-anchor.
            // Decline rather than stamp the parent accession on an unresolved
            // chromosome coordinate, which would be invalid HGVS (#655).
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "cannot re-anchor {} into its NG_/LRG_ parent frame: an endpoint has no \
                     single resolved position (uncertain or compound boundary) (#480/#655)",
                    gv.accession.transcript_accession(),
                ),
            });
        };
        let (Some(a), Some(b)) = (placement.nc_to_parent(nc_lo), placement.nc_to_parent(nc_hi))
        else {
            // Endpoint falls outside the placed genomic span: the parent frame
            // cannot express this coordinate. Decline rather than emit a
            // chromosome coordinate under the parent accession (#655).
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "cannot re-anchor {} into its NG_/LRG_ parent frame: an endpoint falls \
                     outside the placed genomic span (#480/#655)",
                    gv.accession.transcript_accession(),
                ),
            });
        };
        // A minus-strand placement reverses coordinate order.
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let interval = GenomeInterval::new(GenomePos::new(lo), GenomePos::new(hi));
        let edit = if placement.strand == crate::reference::Strand::Minus {
            gv.loc_edit
                .edit
                .map(|inner| transform_edit_for_strand(&inner, crate::reference::Strand::Minus))
        } else {
            gv.loc_edit.edit
        };
        Ok(GenomeVariant {
            accession: gv.accession,
            gene_symbol: gv.gene_symbol,
            loc_edit: LocEdit::with_uncertainty(interval, edit),
        })
    }

    /// If `variant` is a genome variant on an `NG_`/`LRG_` parent whose
    /// chromosomal [`GenomicPlacement`] is known, rewrite it into the chromosome
    /// (`NC_`) frame — the inverse of [`Self::reanchor_genome_to_parent`] (#480).
    ///
    /// cdot indexes and aligns transcripts on `NC_` contigs, so an `NG_`/`LRG_`
    /// genomic *input* (e.g. `NG_012337.1:g.4812_4813insTAC`) cannot be
    /// enumerated or projected in its own frame. Translating it into the `NC_`
    /// frame lets the normal fan-out path find and project onto the overlapping
    /// transcripts.
    ///
    /// Returns the rewritten variant and the original parent accession when a
    /// rewrite happened (so the per-transcript output can be re-framed under the
    /// parent); otherwise returns the input unchanged with `None`.
    fn deanchor_genomic_parent_input(
        &self,
        variant: &HgvsVariant,
    ) -> (HgvsVariant, Option<Accession>) {
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::variant::GenomeVariant;

        // The parent placement is selected by the build the input variant carries
        // (an explicit `NC_*.10`/`.11` or `GRCh37(...)` tag; `None` for a bare
        // `NG_`, which defaults to GRCh38). `None` keeps the prior GRCh38-preferred
        // behavior, so this is inert until an explicit build is present (#653/#713).
        let build = self.build_hint_for_variant(variant);

        // Transcript-coordinate inputs (c./n./r.) carrying an NG_/LRG_
        // genomic_context are taken into the chromosome (NC_) frame here too, so
        // the fan-out below enumerates the OTHER overlapping transcripts
        // (cross-isoform) and frames each result under the parent (#646). Without
        // this they project only onto the input transcript, unframed.
        let cnr_accession = match variant {
            HgvsVariant::Cds(c) => Some(&c.accession),
            HgvsVariant::Tx(t) => Some(&t.accession),
            HgvsVariant::Rna(r) => Some(&r.accession),
            _ => None,
        };
        if let Some(acc) = cnr_accession {
            if let Some(ctx) = acc.genomic_context.as_deref() {
                if let Some(placement) = (ctx.prefix.as_ref() == "NG" || ctx.is_lrg())
                    .then(|| self.provider.genomic_placement_on_build(ctx, build))
                    .flatten()
                {
                    if let Ok(HgvsVariant::Genome(nc_gv)) = self.project_to_genomic_nc(variant) {
                        // `nc_gv` carries the parent (NG_/LRG_) accession but NC_
                        // coordinates. Stamp the placement's own `NC_` accession —
                        // the authoritative chromosome for this parent — so the
                        // projector enumerates transcripts at the locus. This
                        // mirrors the bare-genomic branch below and
                        // `reanchor_genome_to_parent`, which also key off
                        // `placement.nc`; deriving the contig from cdot instead can
                        // disagree when cdot stores a `chr`-style alias rather than
                        // the `NC_` accession.
                        let nc_variant = HgvsVariant::Genome(GenomeVariant {
                            accession: placement.nc.clone(),
                            gene_symbol: None,
                            loc_edit: nc_gv.loc_edit,
                        });
                        return (nc_variant, Some(ctx.clone()));
                    }
                }
            }
        }

        let HgvsVariant::Genome(gv) = variant else {
            return (variant.clone(), None);
        };
        let Some(placement) = self
            .provider
            .genomic_placement_on_build(&gv.accession, build)
        else {
            return (variant.clone(), None);
        };
        let bounds = match (
            gv.loc_edit.location.start.inner(),
            gv.loc_edit.location.end.inner(),
        ) {
            (Some(s), Some(e)) => Some((s.base, e.base)),
            _ => None,
        };
        let Some((p_lo, p_hi)) = bounds else {
            return (variant.clone(), None);
        };
        let (Some(a), Some(b)) = (placement.parent_to_nc(p_lo), placement.parent_to_nc(p_hi))
        else {
            return (variant.clone(), None);
        };
        // A minus-strand placement reverses coordinate order.
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let interval = GenomeInterval::new(GenomePos::new(lo), GenomePos::new(hi));
        // Resolve any range-reference insertion (`ins<start>_<end>`, `delins…`)
        // against the *parent* (`NG_`/`LRG_`) sequence in the parent frame BEFORE
        // de-anchoring, so the de-anchored variant carries a frame-independent
        // literal. The range-ref endpoints are parent-frame coordinates; left
        // unresolved they would be fetched against the `NC_` accession stamped
        // below and come back as `insNNN…` (#652). Best-effort: if the parent
        // sequence can't be fetched, keep the original edit (no regression).
        let parent_resolved_edit = gv.loc_edit.edit.clone().map(|inner| {
            match crate::normalize::rules::canonicalize_insertion_expand(
                &inner,
                // The bare parent accession (strips any genomic_context), matching
                // the accessor `reanchor_genome_to_parent` uses — `full()` would
                // prepend a `CTX(…)` wrapper and fail the provider name lookup.
                &gv.accession.transcript_accession(),
                crate::normalize::rules::InsCoordKind::Direct,
                &self.provider,
            ) {
                // Range reference resolved to a literal — carry the frame-independent bases.
                Ok(Some(resolved)) => resolved,
                // Legitimate no-op: a plain literal / count / non-range payload returns
                // `Ok(None)` (src/normalize/rules.rs); keep the original edit.
                Ok(None) => inner,
                // A genuine resolution failure (out-of-bounds range, non-IUPAC base,
                // provider error). Do NOT silently re-emit the unresolved
                // `ins<start>_<end>` as `insNNN…` (#652) without a trace: surface it so
                // an unresolvable range reference is at least visible, even though the
                // `(HgvsVariant, Option<Accession>)` return type has no `Result` channel
                // to decline through.
                Err(e) => {
                    log::warn!(
                        "deanchor_genomic_parent_input: failed to resolve range-reference \
                         insertion against parent {}: {} — keeping the unresolved edit",
                        gv.accession.transcript_accession(),
                        e,
                    );
                    inner
                }
            }
        });
        let edit = if placement.strand == crate::reference::Strand::Minus {
            parent_resolved_edit
                .map(|inner| transform_edit_for_strand(&inner, crate::reference::Strand::Minus))
        } else {
            parent_resolved_edit
        };
        let parent = gv.accession.clone();
        let nc_variant = HgvsVariant::Genome(GenomeVariant {
            accession: placement.nc.clone(),
            gene_symbol: gv.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(interval, edit),
        });
        (nc_variant, Some(parent))
    }

    /// Project an already-normalized g. variant onto the *curated* set of
    /// overlapping transcripts (#656), skipping re-normalization.
    ///
    /// The result is the curated enumerated set rather than every
    /// cdot-overlapping record: superseded transcript versions are collapsed
    /// (only the highest version per base accession is kept) and predicted
    /// `XM_`/`XR_` models are dropped when a curated `NM_`/`NR_` transcript
    /// covers the same locus; predicted models are kept only when they are the
    /// sole coverage. See [`select_enumerated_transcript_ids`] for the policy.
    ///
    /// Callers that pre-normalize once and then fan-out across transcripts
    /// should use this method.
    pub fn project_normalized_all(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VariantProjection>, FerroError> {
        // Surface a contradictory `--assembly` override like every other public
        // entry point (#715). The internal caller `project_variant_all` already
        // warns before delegating, so it routes through the non-warning
        // `project_normalized_all_inner` to avoid double-warning.
        self.warn_assembly_conflict(variant);
        self.project_normalized_all_inner(variant)
    }

    /// Fan-out projection without the `--assembly`-conflict warning, for callers
    /// (e.g. [`project_variant_all`]) that have already emitted it.
    fn project_normalized_all_inner(
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

        // 3a. Apply the enumeration policy (#656): collapse superseded versions
        //     and prefer curated transcripts over predicted models, so the
        //     enumerated set matches mutalyzer's curated set rather than every
        //     cdot-overlapping record.
        let all_ids: Vec<&str> = projection_result
            .projections
            .iter()
            .map(|p| p.transcript_id.as_str())
            .collect();
        let keep: std::collections::HashSet<&str> = select_enumerated_transcript_ids(&all_ids)
            .into_iter()
            .collect();

        // 4. Project against each kept overlapping transcript (priority order preserved).
        let mut results = Vec::with_capacity(keep.len());
        for tx_proj in &projection_result.projections {
            if !keep.contains(tx_proj.transcript_id.as_str()) {
                continue;
            }
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
            let build_hint = self.build_hint_for_variant(original);
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
                rna: None,
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
            rna: None,
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
    /// Returns `Ok(None)` when no protein consequence applies — most intronic
    /// edits, UTR, a non-coding transcript, or an edit shape with no
    /// in-frame/frameshift prediction. The intronic carve-out: a deletion whose
    /// endpoints are both intronic and bracket one or more complete coding exons
    /// (`whole_exon_deletion_span` returns `Some`) *is* predicted — it routes
    /// through the indel predictor on its clamped exonic CDS span (#498) — so
    /// only non-whole-exon intronic edits return `None`. `Err` is reserved for
    /// genuine failures (unknown reference, etc.), never for "no consequence"
    /// (mirrors the existing contract).
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
        // `p.(Met1?)` on a canonical-`ATG`-initiation transcript, or as `p.?`
        // on a non-AUG-initiation transcript where residue 1 is not `Met`
        // (#771; see the start-codon split at the `affects_init` branch below)
        // (HGVS recommendations/protein/{substitution.md:51,
        // deletion.md:62}; #512). This holds even when the edit's 5' end lies
        // in the 5'UTR (CDS base ≤ 0) — e.g. a boundary-spanning insertion or
        // duplication like `c.-1_2dup` — so the coarse `is_utr` gate must not
        // drop it to "no protein predicted" (#504). Intronic offsets never
        // reach the initiation codon, so an offset disqualifies the edit here.
        let affects_init = !is_intronic
            && cds_start.offset.is_none()
            && cds_end.offset.is_none()
            && affects_initiation_codon(c_edit, cds_start.base, cds_end.base);
        // A deletion that removes one or more complete coding exons (both
        // endpoints intronic, bracketing the exon) has a predictable protein
        // consequence (HGVS protein/deletion.md "one or more exons"; #498).
        // Clamp it to the exonic CDS bases removed; if that yields a span, let
        // it through the intronic gate so the indel predictor runs on it below.
        let whole_exon = whole_exon_deletion_span(c_edit, cds_start, cds_end);
        if (is_intronic && whole_exon.is_none()) || !is_coding || (is_utr && !affects_init) {
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
        // those too (for non-initiation edits) rather than risk wrong-frame
        // translation of the common inconsistent-annotation case. An
        // initiation-codon-affecting edit is the exception — it reports an
        // unknown form (`p.?`) without translating, so it is handled below
        // before this decline (#771).
        //
        // A failure to read the CDS (missing sequence, degenerate coords) is
        // left to the per-edit paths below, which already treat
        // `ProteinSequenceUnavailable` as a non-fatal "no protein" — so we
        // only act on a CDS we could read whose start is not `ATG`. Read just
        // the start codon (not the whole 1–10 kbp CDS) since that is all the
        // guard inspects.
        let tx_for_cds_check =
            self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
        // Read the start codon once; it drives both the initiation-codon-unknown
        // form below and the #625 non-ATG decline. An unreadable CDS (missing
        // sequence / degenerate coords) is left to the per-edit paths — treat it
        // as "valid" here so we neither decline nor mis-pick the `p.?` form.
        let start_is_atg = read_cds_start_codon(&tx_for_cds_check)
            .map(|c| cds_has_valid_start(&c))
            .unwrap_or(true);
        // An initiation-codon-affecting edit short-circuits to an unknown-protein
        // form here, before edit-type dispatch, so boundary-spanning edits whose
        // 5' end is in the 5'UTR (which the per-edit paths reject via their
        // `cds base > 0` guards) are still reported. The in-CDS start-codon cases
        // (e.g. `c.1A>G`) also funnel through here; the per-edit paths keep their
        // own `affects_initiation_codon` guard as defense in depth for direct
        // callers (#504, #512). The consequence is reported WITHOUT translating,
        // so the non-ATG #625 frame concern does not apply: a canonical ATG start
        // gives `p.(Met1?)`; a non-AUG-initiation transcript (#625) gives the
        // plain whole-protein unknown `p.?` (residue 1 is not `Met`, so the
        // `Met1?` form would be wrong) — #771.
        if affects_init {
            let tx_for_codon =
                self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
            let prot = if start_is_atg {
                build_initiator_unknown(&prot_acc, &tx_for_codon)
            } else {
                build_whole_protein_unknown(&prot_acc, &tx_for_codon)
            };
            return Ok(Some(prot));
        }
        // Issue #625 (non-initiation edits only — initiation-codon edits already
        // returned an unknown form above). Decline to translate when the CDS does
        // not begin with `ATG`: `cds_start` is then likely inconsistent with the
        // FASTA and translation would read the wrong frame and fabricate a bogus
        // consequence. This also conservatively declines the rare legitimate
        // non-AUG-initiation transcript (CUG/GUG with a `transl_except`) for any
        // non-initiation edit rather than risk the common inconsistent-annotation
        // case.
        if !start_is_atg {
            log::trace!(
                "declining protein prediction for {transcript_id}: reference CDS does not begin \
                 with an ATG start codon (cds_start is inconsistent with the transcript FASTA)",
            );
            return Ok(None);
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
            // Whole-exon deletion (both endpoints intronic): predict from the
            // clamped exonic CDS span, reusing the indel predictor on a clean
            // exonic deletion of that span (#498).
            NaEdit::Deletion { .. } if whole_exon.is_some() => {
                let (lo, hi) = whole_exon.expect("checked is_some");
                let tx_for_codon =
                    self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
                let ref_bundle =
                    self.cached_ref_translation(cache_variant, transcript_id, &tx_for_codon)?;
                let exonic_del = NaEdit::Deletion { sequence: None, length: None };
                match predict_indel_protein(
                    &tx_for_codon,
                    &ref_bundle,
                    lo,
                    hi,
                    &exonic_del,
                    &prot_acc,
                ) {
                    Ok(pv) => protein = Some(pv),
                    Err(FerroError::UnsupportedProjection { .. })
                    | Err(FerroError::ProteinSequenceUnavailable { .. }) => {}
                    Err(other) => return Err(other),
                }
            }
            NaEdit::Deletion { .. }
            | NaEdit::Insertion { .. }
            | NaEdit::Duplication { .. }
            | NaEdit::Delins { .. }
            | NaEdit::Inversion { .. }
                // Only predict for concrete exonic CDS positions. This arm's
                // own `offset.is_none()` guards exclude intronic offsets — the
                // loosened intronic gate above no longer does, since whole-exon
                // deletions now pass it (they are handled by the arm above).
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

        let rna = self
            .cached_get_transcript_for_variant(normalized, transcript_id)
            .ok()
            .and_then(|tx| crate::project::rna::predict_rna(normalized, &tx));

        Ok(VariantProjection {
            genomic: None,
            coding: Some(normalized.clone()),
            noncoding,
            protein,
            rna,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
        })
    }

    /// Direct n./r.→CDS→p. projection for a bare transcript non-coding (`n.`)
    /// or RNA (`r.`) input (no `genomic_context` parent). The c.→g.→CDS
    /// roundtrip cannot run without a genome alignment, but an `n.` base is a
    /// 1-based transcript position: convert it to a CDS position via the
    /// transcript's `cds_start` (the exon/CIGAR-aware
    /// [`CoordinateMapper::tx_to_cds`]), then run the same protein-consequence
    /// prediction the coding path uses. The resulting projection has
    /// `genomic = None` (no genomic representation is available for a bare-NM_
    /// input) (#506).
    ///
    /// `normalized` must be an [`HgvsVariant::Tx`] or [`HgvsVariant::Rna`].
    fn project_noncoding_direct(
        &self,
        normalized: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        use crate::reference::Strand as RefStrand;

        // Pull the common pieces (accession, gene_symbol, edit, raw start/end
        // as a transcript TxPos) out of the n. / r. input. RNA positions share
        // the same `base`/`offset`/`utr3` shape as transcript positions, so we
        // normalize both to `TxPos`. RNA edits carry `Base::U`; the protein
        // machinery reads the transcript's DNA CDS, so we translate U→T below.
        let (accession, gene_symbol, edit_mu, axis_label, start_tx, end_tx) = match normalized {
            HgvsVariant::Tx(v) => {
                let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "n.", "start")?;
                let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "n.", "end")?;
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    "n.",
                    s,
                    e,
                )
            }
            HgvsVariant::Rna(v) => {
                let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "r.", "start")?;
                let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "r.", "end")?;
                // RnaPos shape mirrors TxPos: base / offset / utr3↔downstream.
                let to_tx = |p: crate::hgvs::location::RnaPos| TxPos {
                    base: p.base,
                    offset: p.offset,
                    downstream: p.utr3,
                };
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    "r.",
                    to_tx(s),
                    to_tx(e),
                )
            }
            _ => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "project_noncoding_direct only accepts n./r. inputs, got {}",
                        normalized.variant_type()
                    ),
                })
            }
        };

        // Mirror the genome path's transcript_id-mismatch guard: an n./r. input
        // must be projected against its own transcript.
        let input_tx = accession.transcript_accession();
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

        // `?` sentinel (base == 0) has no concrete transcript position; reject
        // it up front to match the coding-path contract. (Transcript positions
        // cannot carry chromosome-arm markers — only `CdsPos` has a `special`
        // field — so there is no pter/qter/cen case to handle here.)
        if start_tx.base == 0 || end_tx.base == 0 {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "cannot resolve `?` position sentinel on direct {axis_label}→p. projection"
                ),
            });
        }

        // 3'-downstream (`*N`) positions carry the `*N` *relative* 3'UTR offset
        // in `base` with `downstream = true`. `CoordinateMapper::tx_to_cds`
        // (used below) classifies purely by comparing `base` against the CDS
        // bounds and never reads this flag, so it would misread `n.*5` / `r.*5`
        // as in-transcript position 5 — yielding a bogus `c.` form and a
        // spurious protein prediction. The genome-pivot path reshapes
        // `downstream` into `CdsPos.utr3` (which `cds_to_genome` honors), but the
        // direct path has no exon-aware 3'UTR translation, so reject these up
        // front rather than mis-project them.
        if start_tx.downstream || end_tx.downstream {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "cannot project 3'-downstream (`*N`) positions on direct {axis_label}→p. \
                     projection without a genome alignment"
                ),
            });
        }

        // The transcript record is needed both for its CDS structure (to
        // convert tx→CDS) and for protein metadata. Prefer the cdot view for
        // gene/protein metadata, but the full `Transcript` is what carries the
        // exon table and `cds_start`/`cds_end` used by `tx_to_cds`.
        let tx = self.cached_get_transcript_for_variant(normalized, transcript_id)?;
        // Coding/protein/gene metadata: prefer cdot, fall back to the sequence
        // provider's transcript record. Mirror the coding-path sibling
        // (`project_coding_direct`) by resolving `is_coding` cdot-first too — the
        // two transcript stores can disagree (cdot lists a CDS, the provider
        // record doesn't), and reading `is_coding` from a different source than
        // the sibling would silently drop protein where the sibling predicts it.
        let (is_coding, cdot_protein, gene_symbol_meta) =
            match self.projector.mapper().cdot().get_transcript(transcript_id) {
                Some(t) => (
                    t.cds_start.is_some(),
                    t.protein.clone(),
                    t.gene_name.clone(),
                ),
                None => (
                    tx.cds_start.is_some() && tx.cds_end.is_some(),
                    tx.protein_id.clone(),
                    tx.gene_symbol.clone(),
                ),
            };
        // Carry the input's gene symbol through to the derived axes when it has
        // one, else the metadata gene symbol; mirrors the coding path.
        let gene_symbol = gene_symbol.or(gene_symbol_meta);

        let edit = edit_mu
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::UnsupportedProjection {
                reason: format!("{axis_label} variant has no concrete edit"),
            })?;
        // Translate any RNA `U` bases to DNA `T` so the edit reads against the
        // transcript's DNA CDS. On the transcript sense strand no reverse-
        // complement is needed (the c./n. form already reads in transcript
        // orientation), so we pass `Strand::Plus`, whose only effect is the
        // U→T mapping. This is a no-op for n. (DNA) inputs.
        let c_edit = transform_edit_for_strand(&edit, RefStrand::Plus);

        // Non-coding transcript: there is no CDS to convert into, so there is
        // no protein consequence and no c. form. The n. form is the input
        // itself; report it on the non-coding axis and stop.
        if !is_coding {
            return Ok(VariantProjection {
                genomic: None,
                coding: None,
                noncoding: Some(normalized.clone()),
                protein: None,
                // No coding form on a non-coding transcript; the RNA prediction
                // surface (#485) predicts from the c. form, so it's unavailable here.
                rna: None,
                transcript_id: transcript_id.to_string(),
                gene_symbol,
                is_frameshift: false,
                is_intronic: start_tx.offset.is_some() || end_tx.offset.is_some(),
                is_utr: false,
            });
        }

        // Convert the 1-based transcript positions to CDS positions via the
        // exon/CIGAR-aware mapper (the inverse of `noncoding_from_coding`'s
        // `cds_to_tx`). The intronic offset is carried through.
        let mapper = crate::convert::mapper::CoordinateMapper::new(&tx);
        let cds_start = mapper.tx_to_cds(&start_tx)?;
        let cds_end = mapper.tx_to_cds(&end_tx)?;

        // Flags derived from the resolved CDS positions (no genome): intronic =
        // any offset; UTR = 5'UTR (base ≤ 0) or 3'UTR (`*`). Mirrors the coding
        // path so the protein gate behaves identically.
        let is_intronic = cds_start.offset.is_some() || cds_end.offset.is_some();
        let is_utr = !is_intronic
            && (cds_start.base <= 0 || cds_end.base <= 0 || cds_start.utr3 || cds_end.utr3);

        let protein = self.predict_protein_consequence(
            transcript_id,
            cdot_protein,
            is_coding,
            is_intronic,
            is_utr,
            &c_edit,
            &cds_start,
            &cds_end,
            normalized,
        )?;

        // Build the c. (coding) form from the resolved CDS positions and the
        // (U→T-translated) edit. Best-effort: keep the bare transcript
        // accession so it renders without an NC_* wrapper, matching the
        // genome-pivot `coding` field.
        let coding = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession(transcript_id),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new(CdsInterval::new(cds_start, cds_end), c_edit.clone()),
        });

        let frameshift = is_frameshift(normalized);

        let rna = crate::project::rna::predict_rna(&coding, &tx);
        Ok(VariantProjection {
            genomic: None,
            coding: Some(coding),
            // The non-coding axis is the input itself.
            noncoding: Some(normalized.clone()),
            protein,
            rna,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
        })
    }

    /// Fetch the exon's genome and transcript sequence slices in transcript
    /// 5'→3' orientation, ready to feed [`crate::data::exon_realign`].
    ///
    /// Returns `None` (decline) on any condition where a sequence-aware
    /// correction should not be attempted: the exon already carries a cdot
    /// CIGAR (its alignment is modelled — leave the existing E3007 / arithmetic
    /// behavior untouched), genomic data is unavailable, or a slice cannot be
    /// fetched. The two slices are returned `(genome, tx)` with the genome slice
    /// reverse-complemented on a minus-strand transcript so both read in
    /// transcript orientation.
    fn exon_slices_for_realign(
        &self,
        cdot_tx: &CdotTranscript,
        transcript_id: &str,
        exon_idx: usize,
    ) -> Option<(Vec<u8>, Vec<u8>)> {
        // Only the ungapped-but-divergent class is in scope (#644). Exons cdot
        // models with a CIGAR keep their existing semantics.
        match cdot_tx.exon_cigars.get(exon_idx) {
            Some(None) | None => {}
            Some(Some(_)) => return None,
        }
        if !self.provider.has_genomic_data() {
            return None;
        }
        let exon = cdot_tx.exons.get(exon_idx)?;
        let (genome_start, genome_end, tx_start, tx_end) = (exon[0], exon[1], exon[2], exon[3]);
        if genome_end <= genome_start || tx_end <= tx_start {
            return None;
        }

        // Transcript slice: cdot tx coords are 0-based half-open; the provider's
        // `get_sequence` is also 0-based half-open.
        let tx_seq = self
            .provider
            .get_sequence(transcript_id, tx_start, tx_end)
            .ok()?;
        // Genome slice: cdot genome coords are HGVS-value-based, so HGVS g.X is
        // 0-based X-1. The exon spans HGVS [genome_start, genome_end), i.e.
        // 0-based [genome_start - 1, genome_end - 1).
        let genome_seq = self
            .provider
            .get_genomic_sequence(&cdot_tx.contig, genome_start - 1, genome_end - 1)
            .ok()?;

        let genome_bytes = if cdot_tx.strand == crate::reference::Strand::Minus {
            crate::sequence::reverse_complement(&genome_seq).into_bytes()
        } else {
            genome_seq.into_bytes()
        };
        Some((genome_bytes, tx_seq.into_bytes()))
    }

    /// Apply the sequence-aware exon-indel correction to a genome→CDS mapping
    /// (issue #644).
    ///
    /// `naive` is the `CdsPos` the cdot arithmetic produced for `genome_pos`.
    /// When the containing exon is reported ungapped by cdot but its transcript
    /// and genome sequences differ by an indel, the naive position is off by the
    /// indel size; this re-derives the position from a local realignment. Only
    /// simple exonic positions are corrected — intronic / special / UTR-range
    /// positions (which carry an `offset` or a `special` marker) are returned
    /// unchanged, as is any position the realignment cannot confidently place.
    ///
    /// Returns `AlignmentGap` when the corrected position would fall strictly
    /// inside a genome-only gap discovered by the realignment (a genome base
    /// with no transcript counterpart) — consistent with the #603 E3007
    /// semantics for modelled CIGAR deletions.
    fn correct_cds_for_exon_indel(
        &self,
        cdot_tx: &CdotTranscript,
        transcript_id: &str,
        info: &MappingInfo,
        genome_pos: &GenomePos,
        naive: &CdsPos,
    ) -> Result<CdsPos, FerroError> {
        use crate::data::exon_realign::{correct_genome_offset, ExonOffsetCorrection};

        // Only correct simple exonic positions. Intronic offsets and special
        // sentinels are derived by boundary arithmetic the realignment doesn't
        // model; leave them alone.
        if info.is_intronic || naive.offset.is_some() || naive.special.is_some() {
            return Ok(*naive);
        }
        // The mapping records the 1-based exon number; the parallel exon /
        // exon_cigars index is `number - 1`.
        let exon_number = match info.exon_numbers.first() {
            Some(n) => *n,
            None => return Ok(*naive),
        };
        let exon_idx = (exon_number as usize).saturating_sub(1);
        let exon = match cdot_tx.exons.get(exon_idx) {
            Some(e) => e,
            None => return Ok(*naive),
        };
        let (genome_start, genome_end) = (exon[0], exon[1]);

        let (genome_seq, tx_seq) =
            match self.exon_slices_for_realign(cdot_tx, transcript_id, exon_idx) {
                Some(slices) => slices,
                None => return Ok(*naive),
            };

        // Offset of the queried genome position from the exon's transcript-5'
        // end, measured on the genome axis — the same orientation the slices
        // are in (see `cigar_deletion_gap_at_genome_pos`).
        let genome_offset = match cdot_tx.strand {
            crate::reference::Strand::Plus => {
                if genome_pos.base < genome_start || genome_pos.base >= genome_end {
                    return Ok(*naive);
                }
                (genome_pos.base - genome_start) as usize
            }
            crate::reference::Strand::Minus => {
                if genome_pos.base < genome_start || genome_pos.base >= genome_end {
                    return Ok(*naive);
                }
                (genome_end - 1 - genome_pos.base) as usize
            }
            crate::reference::Strand::Unknown => return Ok(*naive),
        };

        match correct_genome_offset(&genome_seq, &tx_seq, genome_offset) {
            None => Ok(*naive),
            Some(ExonOffsetCorrection::InsideGenomeOnlyGap) => Err(FerroError::AlignmentGap {
                msg: format!(
                    "genomic position {} falls in a transcript-genome alignment gap \
                     (sequence-detected genome-only indel) in exon {}",
                    genome_pos.base, exon_number
                ),
            }),
            Some(ExonOffsetCorrection::Delta(0)) => Ok(*naive),
            Some(ExonOffsetCorrection::Delta(delta)) => {
                // Apply the correction in transcript space and re-derive the
                // CDS position so a correction crossing a CDS/UTR boundary is
                // resolved correctly.
                let naive_tx = (exon[2] as i64) + (genome_offset as i64);
                let corrected_tx = naive_tx + delta;
                if corrected_tx < 0 {
                    return Ok(*naive);
                }
                cdot_tx
                    .cds_pos_from_tx_pos(corrected_tx as u64)
                    .map(Ok)
                    .unwrap_or(Ok(*naive))
            }
        }
    }

    /// Apply the sequence-aware exon-indel correction to a transcript→genome
    /// mapping — the c.→g. mirror of [`Self::correct_cds_for_exon_indel`]
    /// (issue #644).
    ///
    /// `tx_pos` is the 0-based transcript position being mapped and `naive` is
    /// the genome coordinate the cdot arithmetic produced for it. When the
    /// containing exon is reported ungapped by cdot but its sequences differ by
    /// an indel, the naive genome coordinate is off by the indel size; this
    /// re-derives it from the local realignment. Returns `naive` unchanged when
    /// no confident correction applies, and `AlignmentGap` when the transcript
    /// position has no genome counterpart (a transcript-only indel base).
    fn correct_genome_for_exon_indel(
        &self,
        cdot_tx: &CdotTranscript,
        transcript_id: &str,
        tx_pos: u64,
        naive: u64,
    ) -> Result<u64, FerroError> {
        use crate::data::exon_realign::{correct_tx_offset, TxOffsetCorrection};

        let exon = match cdot_tx.exon_for_tx_pos(tx_pos) {
            Some(e) => e,
            None => return Ok(naive),
        };
        let exon_idx = (exon.number as usize).saturating_sub(1);
        let (genome_seq, tx_seq) =
            match self.exon_slices_for_realign(cdot_tx, transcript_id, exon_idx) {
                Some(slices) => slices,
                None => return Ok(naive),
            };

        let tx_offset = (tx_pos - exon.tx_start) as usize;
        match correct_tx_offset(&genome_seq, &tx_seq, tx_offset) {
            None | Some(TxOffsetCorrection::Delta(0)) => Ok(naive),
            Some(TxOffsetCorrection::InsideTxOnlyGap) => Err(FerroError::AlignmentGap {
                msg: format!(
                    "transcript position {} falls in a transcript-genome alignment gap \
                     (sequence-detected transcript-only indel) in exon {}",
                    tx_pos, exon.number
                ),
            }),
            Some(TxOffsetCorrection::Delta(delta)) => {
                // The delta is in the genome axis, measured from the exon's
                // transcript-5' end. Apply it to the naive genome coordinate in
                // the strand's genomic direction (the naive coord already came
                // from `genome_start + tx_offset` on plus, `genome_end - 1 -
                // tx_offset` on minus, so on minus the axis runs the other way).
                let corrected = match cdot_tx.strand {
                    crate::reference::Strand::Plus => naive as i64 + delta,
                    crate::reference::Strand::Minus => naive as i64 - delta,
                    crate::reference::Strand::Unknown => return Ok(naive),
                };
                if corrected < 0 {
                    return Ok(naive);
                }
                Ok(corrected as u64)
            }
        }
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
        // Direct n./r.→CDS→p. path: a bare non-coding (`n.`) or RNA (`r.`)
        // input likewise has no genome alignment. An `n.` base is a 1-based
        // transcript position; convert it to a CDS position via the
        // transcript's `cds_start` and run the same protein prediction the
        // coding path uses (#506).
        match normalized {
            HgvsVariant::Tx(t) if t.accession.genomic_context.is_none() => {
                return self.project_noncoding_direct(normalized, transcript_id);
            }
            HgvsVariant::Rna(r) if r.accession.genomic_context.is_none() => {
                return self.project_noncoding_direct(normalized, transcript_id);
            }
            _ => {}
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
                // Pivot in the chromosome (NC_) frame: genome→CDS/protein
                // derivation below runs through cdot, which is NC-based. The
                // re-anchored (NG_/LRG_) form is the genomic *output* only and
                // would mis-map back through cdot (#480).
                let g = self.project_to_genomic_nc(normalized)?;
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
                Ok(res) => {
                    // Sequence-aware correction (#644): when cdot reports the
                    // containing exon ungapped but the transcript and genome
                    // FASTAs genuinely differ by an indel, the naive arithmetic
                    // `res.variant` is off by the indel's size. Re-derive the
                    // CDS position from a local realignment of the exon's two
                    // sequences. A no-op (and unchanged behavior) when the
                    // sequences agree, the exon already carries a CIGAR, or the
                    // alignment is not confident enough to correct.
                    let corrected = self.correct_cds_for_exon_indel(
                        cdot_tx,
                        transcript_id,
                        &res.info,
                        gp,
                        &res.variant,
                    )?;
                    Ok((corrected, res.info))
                }
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

        // `.genomic` is the canonical g. *output*. For a g. input it is the
        // (normalized) input itself, emitted as-is. For a c./n./r. input it is
        // the chromosome-frame pivot re-anchored into an `NG_`/`LRG_` parent's
        // own frame (#480) — the derivation above already ran on the `NC_` pivot
        // (`genome_variant`), which must stay in the chromosome frame for cdot.
        // When the parent frame cannot be resolved (no placement, an endpoint
        // outside the placed span, or an uncertain endpoint) the genomic axis is
        // dropped to `None` rather than stamping the parent accession on a
        // chromosome coordinate (invalid HGVS): unlike the g.-only
        // `project_to_genomic`, this path returns a multi-axis projection, so
        // dropping the unframable genomic axis preserves the still-valid
        // coding/noncoding/protein axes (#655/#702, mirroring `frame_projection_owned`).
        let genomic = match normalized {
            HgvsVariant::Genome(_) => Some(projected_genome),
            _ => match projected_genome {
                HgvsVariant::Genome(gv) => self
                    .reanchor_genome_output(gv, self.build_hint_for_variant(normalized))
                    .ok()
                    // #737: re-normalize the genomic *output* in its own
                    // (parent/`NC_`) frame. The pivot was 3'-normalized in
                    // transcript space; for a minus-strand transcript that is the
                    // genome's 5' end, so without a genomic-frame renormalization
                    // a shiftable del/dup/ins lands at the wrong (5') end on the
                    // genome. Normalizing the re-anchored variant against the
                    // output accession's own sequence yields the spec-canonical
                    // genomic-3' form. Falls back to the un-normalized re-anchored
                    // form if normalization fails (e.g. no reference bases); the
                    // dropped error is `trace!`-logged so a regression is
                    // observable rather than silently masked (the fallback
                    // swallows *any* normalize error, not only missing bases).
                    .map(|reanchored| {
                        self.normalizer.normalize(&reanchored).unwrap_or_else(|e| {
                            log::trace!(
                                "genomic-frame renormalization failed for {reanchored}; \
                                 emitting un-normalized re-anchored form: {e}"
                            );
                            reanchored
                        })
                    }),
                other => Some(other),
            },
        };

        // Predicted RNA consequence (#485): derived best-effort from the c. form
        // and the transcript sequence; `None` when not representable.
        let rna = self
            .cached_get_transcript_for_variant(&cache_variant, transcript_id)
            .ok()
            .and_then(|tx| crate::project::rna::predict_rna(&coding, &tx));

        Ok(VariantProjection {
            genomic,
            coding: Some(coding),
            noncoding,
            protein,
            rna,
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

/// Split a transcript accession into its version-stripped base and version:
/// `"NM_000532.5" → ("NM_000532", Some(5))`, `"NM_000532" → ("NM_000532", None)`.
fn split_accession_version(id: &str) -> (&str, Option<u32>) {
    match id.rsplit_once('.') {
        Some((base, ver)) if !ver.is_empty() && ver.bytes().all(|b| b.is_ascii_digit()) => {
            (base, ver.parse::<u32>().ok())
        }
        _ => (id, None),
    }
}

/// Whether a transcript id is a *predicted* RefSeq model (`XM_`/`XR_`), as
/// opposed to a curated `NM_`/`NR_` transcript.
fn is_predicted_model(id: &str) -> bool {
    id.starts_with("XM_") || id.starts_with("XR_")
}

/// Apply the `project_*_all` enumeration policy (#656) to a priority-ordered
/// list of overlapping transcript ids:
///
/// 1. **Collapse superseded versions** — keep only the highest version per base
///    accession (drop `NM_000532.4` when `NM_000532.5` overlaps the same locus).
/// 2. **Prefer curated transcripts** — drop predicted `XM_`/`XR_` models when a
///    curated `NM_`/`NR_` transcript also covers the locus (matching mutalyzer's
///    curated set), but keep predicted models when they are the *sole* coverage,
///    so a locus a predicted model would have covered never returns empty.
///
/// Input order (clinical priority) is preserved among the kept ids.
pub(crate) fn select_enumerated_transcript_ids<'a>(ids: &[&'a str]) -> Vec<&'a str> {
    use std::collections::HashMap;
    // 1. Highest version seen per base accession.
    let mut max_ver: HashMap<&str, u32> = HashMap::new();
    for id in ids {
        let (base, ver) = split_accession_version(id);
        if let Some(v) = ver {
            max_ver
                .entry(base)
                .and_modify(|m| *m = (*m).max(v))
                .or_insert(v);
        }
    }
    let highest: Vec<&str> = ids
        .iter()
        .copied()
        .filter(|id| {
            let (base, ver) = split_accession_version(id);
            match ver {
                Some(v) => max_ver.get(base) == Some(&v),
                // A bare (unversioned) accession is superseded by any versioned
                // form of the same base; keep it only when no versioned id for
                // that base was seen.
                None => !max_ver.contains_key(base),
            }
        })
        .collect();
    // 2. Curated-preferred, never-empty.
    if highest.iter().any(|id| !is_predicted_model(id)) {
        highest
            .into_iter()
            .filter(|id| !is_predicted_model(id))
            .collect()
    } else {
        highest
    }
}

#[cfg(test)]
mod enumeration_policy_tests {
    use super::select_enumerated_transcript_ids;

    #[test]
    fn collapses_superseded_versions() {
        // The #656 example: both versions of two transcripts → keep the highest.
        let ids = [
            "NM_000532.4",
            "NM_000532.5",
            "NM_001178014.1",
            "NM_001178014.2",
        ];
        assert_eq!(
            select_enumerated_transcript_ids(&ids),
            vec!["NM_000532.5", "NM_001178014.2"]
        );
    }

    #[test]
    fn drops_predicted_models_when_curated_present() {
        let ids = ["NM_000532.5", "XM_011512873.2", "XM_005247508.1"];
        assert_eq!(select_enumerated_transcript_ids(&ids), vec!["NM_000532.5"]);
    }

    #[test]
    fn keeps_predicted_models_when_sole_coverage() {
        // No curated transcript overlaps → keep predicted models (never empty).
        let ids = ["XM_005247508.1", "XM_011512873.2"];
        assert_eq!(
            select_enumerated_transcript_ids(&ids),
            vec!["XM_005247508.1", "XM_011512873.2"]
        );
    }

    #[test]
    fn collapses_predicted_versions_then_drops_when_curated_present() {
        let ids = ["NM_000532.5", "XM_011512873.1", "XM_011512873.2"];
        assert_eq!(select_enumerated_transcript_ids(&ids), vec!["NM_000532.5"]);
    }

    #[test]
    fn preserves_priority_order_among_kept() {
        let ids = ["NM_B.1", "NM_A.2", "NM_A.1"];
        assert_eq!(
            select_enumerated_transcript_ids(&ids),
            vec!["NM_B.1", "NM_A.2"]
        );
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
        // A second transcript with a NON-ATG start codon ("CTG…"), for the #771
        // initiation-codon-on-a-non-AUG-transcript path. Bare-NM_ direct
        // projection needs only the provider transcript + cdot version check.
        cdot.add_transcript(
            "NM_NOATG.1".to_string(),
            CdotTranscript {
                gene_name: Some("NOATGGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[2000, 2009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_NOATG.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        // Non-ATG transcript: "CTGCGCTAA" (CTG start), cds_start=1.
        provider.add_transcript(Transcript {
            id: "NM_NOATG.1".to_string(),
            gene_symbol: Some("NOATGGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("CTGCGCTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(2000),
            genomic_end: Some(2008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
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

    /// Like [`make_two_transcript_setup`] but seeds three transcripts on the
    /// *same* locus to exercise the #656 enumeration policy end-to-end through
    /// the `project_*_all` fan-out:
    /// - `NM_TX1.1` and `NM_TX1.2` — superseded/current versions of one curated
    ///   base accession (only `.2` should survive the version collapse);
    /// - `XM_TX9.1` — a predicted model that should be dropped because a curated
    ///   transcript covers the same locus.
    ///
    /// So a curated set of exactly `{NM_TX1.2}` is expected, versus the three
    /// overlapping records cdot reports.
    #[cfg(test)]
    fn make_curated_enumeration_setup() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        for id in ["NM_TX1.1", "NM_TX1.2", "XM_TX9.1"] {
            cdot.add_transcript(
                id.to_string(),
                CdotTranscript {
                    gene_name: Some("GENE1".to_string()),
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
        }
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        for id in ["NM_TX1.1", "NM_TX1.2", "XM_TX9.1"] {
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

    /// Like [`make_two_transcript_setup`] but with the transcripts placed on an
    /// `NC_` accession contig (`NC_000001.11`) rather than the `chr1` alias, so
    /// an `NG_` `GenomicPlacement` (whose `nc` field is an `NC_` accession) maps
    /// onto the same contig the projector indexes. Used to exercise the
    /// NG_-genomic-input fan-out path (#480 inverse).
    #[cfg(test)]
    fn make_nc_two_transcript_setup() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        for (id, np) in [("NM_TX1.1", "NP_TX1.1"), ("NM_TX2.1", "NP_TX2.1")] {
            cdot.add_transcript(
                id.to_string(),
                CdotTranscript {
                    gene_name: Some("GENE1".to_string()),
                    contig: "NC_000001.11".to_string(),
                    strand: ProvStrand::Plus,
                    exons: vec![[1000, 1009, 0, 9]],
                    cds_start: Some(0),
                    cds_end: Some(9),
                    gene_id: None,
                    protein: Some(np.to_string()),
                    exon_cigars: Vec::new(),
                },
            );
        }
        let projector = Projector::new(cdot).with_mane(vec!["NM_TX2.1".to_string()], vec![]);

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
                chromosome: Some("NC_000001.11".to_string()),
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
        provider.add_genomic_sequence(
            "NC_000001.11",
            format!("{}{}{}", prefix, "ATGCGCTAA", suffix),
        );
        // The NG_ parent's own sequence (NG_900.1 base 1..9 == NC 1000..1008),
        // so normalization of an NG_-relative input resolves its reference.
        provider.add_genomic_sequence("NG_900.1", "ATGCGCTAA".to_string());
        (projector, provider)
    }

    /// A bare `NG_<n>:g.` input must enumerate the transcripts overlapping that
    /// NG_ position by translating the NG_-relative coordinate into the
    /// chromosome frame via the parent's [`GenomicPlacement`] (#480 inverse),
    /// then projecting onto each overlapping transcript. Pre-fix
    /// `extract_contig_and_pos` used the bare NG_ accession as the contig, which
    /// the projector cannot map, so `project_variant_all` returned `[]`.
    #[test]
    fn project_variant_all_enumerates_transcripts_for_ng_genomic_input() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        // NG_900.1 placed on NC_000001.11, plus strand, NG_ base 1 == nc 1000.
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        // NG_900.1:g.4C>A — NG_ base 4 == nc 1003 == c.4 of the CDS.
        let results = vp
            .project_all("NG_900.1:g.4C>A")
            .expect("project_all should enumerate transcripts for an NG_ genomic input");

        assert_eq!(
            results.len(),
            2,
            "expected both transcripts overlapping the NG_ position, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// #637: a legacy gene-model selector on a genomic reference
    /// (`NG_900.1(GENE1_v001):c.…`) resolves to the gene's transcript via the
    /// normalize-first step in projection, so it no longer fails with
    /// "Reference not found: NG_900.1" — it projects onto the resolved `NM_`.
    #[test]
    fn legacy_gene_model_selector_projects_via_normalize_rewrite() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        provider.add_legacy_gene_model("GENE1", "NM_TX2.1");
        let vp = VariantProjector::new(projector, provider);

        let proj = vp
            .project("NG_900.1(GENE1_v001):c.4C>A", "NM_TX2.1")
            .expect("legacy gene-model selector must resolve and project, not error on NG_");
        let coding = proj
            .coding
            .as_ref()
            .expect("coding axis present")
            .to_string();
        assert!(
            coding.contains("NM_TX2.1"),
            "the GENE1_v001 selector must resolve to the transcript: {coding}"
        );
    }

    /// #652: a range-reference insertion (`ins<start>_<end>`) on an `NG_`
    /// genomic input must resolve to the literal parent bases — the range-ref
    /// endpoints are in the `NG_` parent frame and, if left unresolved through
    /// the de-anchor into the `NC_` frame, would be fetched against the `NC_`
    /// accession and come back as `insNNN…`.
    #[test]
    fn deanchor_resolves_ng_range_reference_insertion_to_literal() {
        let mut provider = MockProvider::new();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 30000,
                strand: crate::reference::Strand::Plus,
            },
        );
        // NG_900.1 parent sequence: bases 4300..=4320 (1-based) are the insert.
        let insert = "GTCCTGTGCTCATTATCTGGC"; // 21 bases
        let mut seq = vec![b'A'; 4299];
        seq.extend_from_slice(insert.as_bytes());
        provider.add_genomic_sequence("NG_900.1", String::from_utf8(seq).unwrap());

        let vp = VariantProjector::new(Projector::new(CdotMapper::new()), provider);
        let variant = crate::parse_hgvs("NG_900.1:g.5207_5208ins4300_4320").unwrap();
        let (deanchored, parent) = vp.deanchor_genomic_parent_input(&variant);

        assert_eq!(
            parent.as_ref().map(|a| a.full()),
            Some("NG_900.1".to_string())
        );
        let rendered = deanchored.to_string();
        assert!(
            rendered.contains(&format!("ins{insert}")),
            "expected the literal parent bases, got: {rendered}"
        );
        assert!(
            !rendered.contains("insN"),
            "must not emit insNNN (unresolved range reference): {rendered}"
        );
    }

    /// #652, minus-strand placement: the range reference is resolved to the
    /// literal parent bases in the parent frame, then the existing strand
    /// transform reverse-complements that literal for the `NC_` frame (the same
    /// path a literal `insACGT` already takes), never `insNNN`.
    #[test]
    fn deanchor_resolves_ng_range_reference_insertion_minus_strand() {
        let mut provider = MockProvider::new();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 30000,
                strand: crate::reference::Strand::Minus,
            },
        );
        let insert = "GTCCTGTGCTCATTATCTGGC"; // parent bases 4300..=4320
        let revcomp = "GCCAGATAATGAGCACAGGAC"; // reverse-complement of `insert`
        let mut seq = vec![b'A'; 4299];
        seq.extend_from_slice(insert.as_bytes());
        provider.add_genomic_sequence("NG_900.1", String::from_utf8(seq).unwrap());

        let vp = VariantProjector::new(Projector::new(CdotMapper::new()), provider);
        let variant = crate::parse_hgvs("NG_900.1:g.5207_5208ins4300_4320").unwrap();
        let (deanchored, _) = vp.deanchor_genomic_parent_input(&variant);
        let rendered = deanchored.to_string();
        assert!(
            rendered.contains(&format!("ins{revcomp}")),
            "expected the reverse-complemented literal, got: {rendered}"
        );
        assert!(
            !rendered.contains("insN"),
            "must not emit insNNN: {rendered}"
        );
    }

    /// #652, `delins` range reference: `canonicalize_insertion_expand` and the
    /// production doc-comment explicitly handle `delins<start>_<end>`, so a
    /// `delins`-shaped range-ref input must also resolve to the literal parent
    /// bases (never `insNNN…`), exercising the named-but-otherwise-untested edit
    /// kind alongside the plain `Insertion` cases above.
    #[test]
    fn deanchor_resolves_ng_range_reference_delins_to_literal() {
        let mut provider = MockProvider::new();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 30000,
                strand: crate::reference::Strand::Plus,
            },
        );
        // NG_900.1 parent sequence: bases 4300..=4320 (1-based) are the insert.
        let insert = "GTCCTGTGCTCATTATCTGGC"; // 21 bases
        let mut seq = vec![b'A'; 4299];
        seq.extend_from_slice(insert.as_bytes());
        provider.add_genomic_sequence("NG_900.1", String::from_utf8(seq).unwrap());

        let vp = VariantProjector::new(Projector::new(CdotMapper::new()), provider);
        let variant = crate::parse_hgvs("NG_900.1:g.5207_5208delins4300_4320").unwrap();
        let (deanchored, parent) = vp.deanchor_genomic_parent_input(&variant);

        assert_eq!(
            parent.as_ref().map(|a| a.full()),
            Some("NG_900.1".to_string())
        );
        let rendered = deanchored.to_string();
        assert!(
            rendered.contains(&format!("delins{insert}")),
            "expected the literal parent bases, got: {rendered}"
        );
        assert!(
            !rendered.contains("insN"),
            "must not emit insNNN (unresolved range reference): {rendered}"
        );
    }

    /// For an `NG_`-genomic input, each per-transcript coding and protein
    /// description must be framed under the original `NG_` parent
    /// (`NG_900.1(NM_TX1.1):c.…` / `NG_900.1(NP_TX1.1):p.…`), matching the
    /// mutalyzer corpus, rather than the bare transcript or the de-anchored
    /// `NC_` accession used internally for projection.
    #[test]
    fn project_variant_all_frames_coding_protein_under_ng_parent() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("NG_900.1:g.4C>A")
            .expect("project_all should succeed for an NG_ genomic input");

        let coding_strings: Vec<String> = results
            .iter()
            .filter_map(|r| r.coding.as_ref().map(|c| c.to_string()))
            .collect();
        assert!(
            !coding_strings.is_empty(),
            "expected at least one coding description"
        );
        for c in &coding_strings {
            assert!(
                c.starts_with("NG_900.1("),
                "coding description should be framed under the NG_ parent, got {c:?}"
            );
        }
        for r in &results {
            if let Some(p) = r.protein.as_ref().map(|p| p.to_string()) {
                assert!(
                    p.starts_with("NG_900.1("),
                    "protein description should be framed under the NG_ parent, got {p:?}"
                );
            }
            // The genomic description is re-anchored back into the NG_ frame
            // (NG_ base 4 == nc 1003), so it reads `NG_900.1:g.4…`.
            if let Some(g) = r.genomic.as_ref().map(|g| g.to_string()) {
                assert!(
                    g.starts_with("NG_900.1:g.4"),
                    "genomic description should be in the NG_ frame, got {g:?}"
                );
            }
        }
    }

    /// A transcript-coordinate input carrying an `NG_` `genomic_context`
    /// (`NG_900.1(NM_TX1.1):c.4C>A`) must, like a bare `NG_:g.` input, enumerate
    /// the overlapping transcripts and frame each description under the `NG_`
    /// parent — not return only the input transcript, bare (#646 / Mode A).
    #[test]
    fn project_variant_all_frames_coding_input_with_ng_context() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::CdsVariant;
        // NG_900.1(NM_TX1.1):c.4C>A — c. input carrying an NG_ genomic_context.
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NG",
                "900",
                Some(1),
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
        let v = HgvsVariant::Cds(cds);
        let results = vp
            .project_variant_all(&v)
            .expect("project_all should succeed for a c. input with NG_ context");

        // Both overlapping transcripts are enumerated.
        let txs: Vec<&str> = results.iter().map(|r| r.transcript_id.as_str()).collect();
        assert!(
            txs.contains(&"NM_TX1.1") && txs.contains(&"NM_TX2.1"),
            "expected both transcripts enumerated, got {txs:?}"
        );
        // Every coding description is framed under the NG_ parent.
        for r in &results {
            if let Some(c) = r.coding.as_ref().map(|c| c.to_string()) {
                assert!(
                    c.starts_with("NG_900.1("),
                    "coding should be framed under the NG_ parent, got {c:?}"
                );
            }
        }
    }

    /// Coexistence of the #646 cross-isoform enumeration with the #655 graceful
    /// decline. A `c.` input carrying an `NG_` `genomic_context` is de-anchored
    /// into the chromosome frame and enumerated against the overlapping
    /// transcripts, each re-framed under the parent (#646). Here the registered
    /// placement spans a region disjoint from where the transcripts actually map
    /// (chr 1003), so the genomic axis cannot be re-anchored into the parent
    /// frame. Rather than failing the whole projection (the #646 enumerate path)
    /// or stamping the parent accession on an out-of-frame chromosome coordinate
    /// (invalid HGVS, #655), `project_variant_all` returns the transcript-relative
    /// coding/protein forms framed under the parent and drops only the genomic
    /// axis (`None`). This locks in that #646 (enumerate) and #655 (graceful
    /// decline) coexist: `frame_projection_owned` degrades the unframable axis
    /// via `reanchor_genome_to_parent` instead of hard-erroring.
    #[test]
    fn project_variant_all_drops_genomic_axis_but_still_enumerates_under_ng_parent() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        // Placement exists (so de-anchor enumerates) but its span [2000, 2010]
        // is disjoint from chr 1003 where the transcripts map, so the genomic
        // re-anchor `nc_to_parent(1003)` declines and the axis is dropped.
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 2000,
                nc_end: 2010,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::CdsVariant;
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NG",
                "900",
                Some(1),
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
        let v = HgvsVariant::Cds(cds);
        let results = vp
            .project_variant_all(&v)
            .expect("enumerate path must degrade the unframable genomic axis, not hard-error");

        // Enumeration still finds the overlapping transcripts (#646 path alive).
        let txs: Vec<&str> = results.iter().map(|r| r.transcript_id.as_str()).collect();
        assert!(
            txs.contains(&"NM_TX1.1") && txs.contains(&"NM_TX2.1"),
            "expected both transcripts enumerated, got {txs:?}"
        );
        for r in &results {
            // The genomic axis is dropped — it cannot be expressed in the parent
            // frame — rather than emitting a chromosome coordinate under NG_900.1.
            assert!(
                r.genomic.is_none(),
                "genomic axis should be dropped when it cannot be re-anchored, got: {:?}",
                r.genomic
            );
            // The coding axis survives, re-framed under the NG_ parent.
            if let Some(c) = r.coding.as_ref().map(|c| c.to_string()) {
                assert!(
                    c.starts_with("NG_900.1("),
                    "coding should be framed under the NG_ parent, got {c:?}"
                );
            }
        }
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

    // ------------------------------------------------------------------
    // Sequence-aware exon-indel correction (issue #644)
    // ------------------------------------------------------------------

    /// Build a projector + provider for a plus-strand transcript whose cdot
    /// record is **ungapped** (`exon_cigars` empty) but whose transcript and
    /// genome sequences genuinely differ by an indel — the cdot 0.2.32 /
    /// `NM_000532.5` (PCCB) failure shape, reduced to a CI-runnable fixture.
    ///
    /// Single exon, genome HGVS `[1000, 1012)` (12 bp) ↔ tx `[0, 11)` (11 bp),
    /// `cds_start = 0`. The true alignment (recoverable from the sequences) is:
    /// 2 genome-only bases at the exon start, then 10 matched bases, then 1
    /// transcript-only base at the end (net genome − tx = +1). So a genomic
    /// position in the matched region maps to a transcript position 2 lower than
    /// the naive cdot arithmetic would give.
    fn make_ungapped_indel_provider_and_projector() -> (Projector, MockProvider) {
        let common = "ATGCGCTAAC"; // 10 bases, matched region
        let genome_exon = format!("CA{common}"); // 12 bp: 2 genome-only + 10 matched
        let tx_exon = format!("{common}T"); // 11 bp: 10 matched + 1 tx-only
        debug_assert_eq!(genome_exon.len(), 12);
        debug_assert_eq!(tx_exon.len(), 11);

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_GAP644.1".to_string(),
            CdotTranscript {
                gene_name: Some("GAPGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // genome [1000, 1012) (HGVS-value), tx [0, 11) — lengths match
                // the FASTAs, but cdot records NO CIGAR (the bug shape).
                exons: vec![[1000, 1012, 0, 11]],
                cds_start: Some(0),
                cds_end: Some(11),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_GAP644.1".to_string(),
            gene_symbol: Some("GAPGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some(tx_exon.clone()),
            cds_start: Some(1), // 1-based per Transcript convention
            cds_end: Some(11),
            exons: vec![Exon::new(1, 1, 11)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1011),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Genomic contig: 999 N's (so HGVS g.1000 == 0-based index 999) + the
        // exon sequence + trailing N's.
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(50);
        provider.add_genomic_sequence("chr1", format!("{prefix}{genome_exon}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_corrects_genome_to_cds_across_unmodelled_indel() {
        // The reproducer class: a g. SNV in the matched region must project to a
        // c. coordinate corrected by the 2-bp genome-only prefix.
        //
        // First matched genome base is HGVS g.1002 (0-based exon offset 2). It
        // maps to tx offset 0 → c.1. Without the correction the naive arithmetic
        // would give tx offset 2 → c.3.
        let (projector, provider) = make_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // g.1006 is exon offset 6 → naive tx 6 (c.7); corrected tx 4 (c.5).
        let result = vp
            .project("chr1:g.1006A>T", "NM_GAP644.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").to_string();
        assert!(
            c.contains(":c.5"),
            "expected corrected ':c.5' (not naive ':c.7') in '{c}'"
        );
    }

    #[test]
    fn project_first_matched_base_maps_to_c1() {
        let (projector, provider) = make_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1002 is the first matched base → c.1.
        let result = vp
            .project("chr1:g.1002A>T", "NM_GAP644.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").to_string();
        assert!(c.contains(":c.1A>T"), "expected ':c.1A>T' in '{c}'");
    }

    /// Non-coding (`cds_start`/`cds_end == None`) sibling of
    /// `make_ungapped_indel_provider_and_projector`: same 2-bp genome-only
    /// prefix + 10-bp matched + 1-bp tx-only exon shape, but the transcript
    /// carries no CDS so it projects on the `n.` axis.
    fn make_noncoding_ungapped_indel_provider_and_projector() -> (Projector, MockProvider) {
        let common = "ATGCGCTAAC"; // 10 bases, matched region
        let genome_exon = format!("CA{common}"); // 12 bp: 2 genome-only + 10 matched
        let tx_exon = format!("{common}T"); // 11 bp: 10 matched + 1 tx-only

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NR_GAP644.1".to_string(),
            CdotTranscript {
                gene_name: Some("GAPGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1012, 0, 11]],
                // Non-coding: no CDS bounds.
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
            id: "NR_GAP644.1".to_string(),
            gene_symbol: Some("GAPGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some(tx_exon.clone()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 11)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1011),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(50);
        provider.add_genomic_sequence("chr1", format!("{prefix}{genome_exon}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_to_genomic_noncoding_corrects_across_unmodelled_indel() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{Accession, TxVariant};

        // The outbound n.→g. path must apply the same #644 sequence-aware
        // correction the inbound g.→n. path (`correct_cds_for_exon_indel`, which
        // is coding-agnostic) applies, so the two round-trip. Without the
        // correction in the non-coding `else` branch of `map_pos`, n.5 would map
        // to the naive g.1004; with it, the 2-bp genome-only prefix shifts it to
        // g.1006 (tx offset 4 → tx-oriented genome offset 6 → HGVS g.1006).
        let (projector, provider) = make_noncoding_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        let n = TxVariant {
            accession: parse_accession("NR_GAP644.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GAPGENE".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(5)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::T,
                },
            ),
        };
        let g = vp
            .project_to_genomic(&HgvsVariant::Tx(n))
            .expect("n. → g. projection should succeed");
        assert!(
            g.to_string().contains("g.1006"),
            "expected corrected 'g.1006' (not naive 'g.1004'), got '{g}'"
        );
    }

    /// Minus-strand sibling of `make_ungapped_indel_provider_and_projector`.
    /// The transcript-oriented exon is the same shape (2 genome-only prefix +
    /// 10 matched + 1 tx-only), but the transcript reads off the minus strand,
    /// so the genome FASTA stores the reverse complement of the
    /// transcript-oriented exon. This is the only fixture that drives the
    /// minus-strand axis-flip arithmetic in `correct_cds_for_exon_indel`
    /// (`genome_end - 1 - genome_pos.base`) and `correct_genome_for_exon_indel`
    /// (`naive - delta`).
    fn make_minus_ungapped_indel_provider_and_projector() -> (Projector, MockProvider) {
        // Transcript-oriented exon (5'→3' as the transcript reads it):
        let common = "ATGCGCTAAC"; // 10 matched bases
        let tx_oriented_genome_exon = format!("CA{common}"); // 12 bp: 2 genome-only + 10 matched
        let tx_exon = format!("{common}T"); // 11 bp: 10 matched + 1 tx-only
                                            // The genome FASTA stores the minus-strand (reference) orientation:
        let genome_fasta_exon = crate::sequence::reverse_complement(&tx_oriented_genome_exon);

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_GAP644M.1".to_string(),
            CdotTranscript {
                gene_name: Some("GAPGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Minus,
                // genome HGVS [1000, 1012); tx [0, 11). No CIGAR (the bug shape).
                exons: vec![[1000, 1012, 0, 11]],
                cds_start: Some(0),
                cds_end: Some(11),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_GAP644M.1".to_string(),
            gene_symbol: Some("GAPGENE".to_string()),
            strand: TxStrand::Minus,
            sequence: Some(tx_exon.clone()),
            cds_start: Some(1),
            cds_end: Some(11),
            exons: vec![Exon::new(1, 1, 11)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1011),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(50);
        provider.add_genomic_sequence("chr1", format!("{prefix}{genome_fasta_exon}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_minus_strand_corrects_and_roundtrips_across_unmodelled_indel() {
        use crate::hgvs::variant::Accession;

        // On minus strand the transcript-5' end is the exon's genome-3' end
        // (HGVS g.1011). The first matched base is tx-oriented genome offset 2
        // → HGVS g.1009 → c.1. A SNV at HGVS g.1005 is tx-oriented genome
        // offset 6 → naive tx 6 (c.7); the 2-bp genome-only prefix corrects it
        // to tx 4 (c.5).
        let (projector, provider) = make_minus_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // g.→c. correction (drives the minus-strand `correct_cds_for_exon_indel`
        // axis flip). Base at g.1005 is 'C' in reference orientation.
        let result = vp
            .project("chr1:g.1005C>A", "NM_GAP644M.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").clone();
        assert!(
            c.to_string().contains(":c.5"),
            "expected corrected ':c.5' (not naive ':c.7') in '{c}'"
        );

        // c.→g. round-trip (drives the minus-strand
        // `correct_genome_for_exon_indel` `naive - delta` branch).
        let c_with_parent = match c {
            HgvsVariant::Cds(mut cds) => {
                cds.accession =
                    cds.accession
                        .with_genomic_context(Accession::new("NC", "000001", Some(11)));
                HgvsVariant::Cds(cds)
            }
            other => panic!("expected a Cds (c.) variant, got {other:?}"),
        };
        let g = vp
            .project_to_genomic(&c_with_parent)
            .expect("c. → g. round-trip should succeed");
        assert!(
            g.to_string().contains("g.1005"),
            "expected round-trip back to 'g.1005', got '{g}'"
        );
    }

    #[test]
    fn project_inside_genome_only_gap_errors_alignment_gap() {
        // A variant on a genome-only base (the 2-bp prefix) has no transcript
        // counterpart and must surface as AlignmentGap, not a wrong coordinate.
        let (projector, provider) = make_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1000 / g.1001 are the genome-only prefix bases.
        let err = vp
            .project("chr1:g.1000C>A", "NM_GAP644.1")
            .expect_err("expected AlignmentGap for a genome-only base");
        assert!(
            matches!(err, FerroError::AlignmentGap { .. }),
            "expected AlignmentGap, got {err:?}"
        );
    }

    #[test]
    fn project_identical_sequences_unaffected_by_correction() {
        // Control: the normal fixture (transcript == genome sequence) must be
        // byte-identical to its pre-#644 behavior — no spurious correction.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "NM_TEST.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").to_string();
        assert!(
            c.contains(":c.4C>A"),
            "expected uncorrected ':c.4C>A' in '{c}'"
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
        // cdot genome coords are HGVS-value-based (HGVS g.X == 0-based index
        // X-1), matching every other fixture here, so the exon's genome_start
        // 1000 lands the CDS at 0-based index 999 (999 N's of prefix). The
        // sequence-aware exon-indel correction (#644) reads the genome FASTA
        // through these cdot coords, so the FASTA must be placed consistently or
        // the realign would fabricate a phantom 1-bp indel against the identical
        // transcript sequence.
        let prefix = "N".repeat(999);
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

    #[test]
    fn project_all_applies_curated_enumeration_policy() {
        // End-to-end: the #656 curation must be wired into the fan-out loop, not
        // just unit-tested on `select_enumerated_transcript_ids`. cdot reports
        // three overlapping records (NM_TX1.1, NM_TX1.2, XM_TX9.1) but the
        // curated set is exactly {NM_TX1.2}: the superseded .1 is collapsed and
        // the predicted XM_ model is dropped in favor of the curated transcript.
        let (projector, provider) = make_curated_enumeration_setup();
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("chr1:g.1003C>A")
            .expect("project_all should succeed");
        let ids: Vec<&str> = results.iter().map(|p| p.transcript_id.as_str()).collect();
        assert_eq!(
            ids,
            vec!["NM_TX1.2"],
            "expected only the curated, current-version transcript, got {:?}",
            ids
        );
    }

    #[test]
    fn project_normalized_all_applies_curated_enumeration_policy() {
        // Same curation, driven through the pre-normalized fan-out entrypoint.
        let (projector, provider) = make_curated_enumeration_setup();
        let vp = VariantProjector::new(projector, provider);

        let variant = crate::parse_hgvs("chr1:g.1003C>A").expect("parse should succeed");
        let normalized = vp.normalizer.normalize(&variant).expect("normalize failed");
        let results = vp
            .project_normalized_all(&normalized)
            .expect("project_normalized_all should succeed");
        let ids: Vec<&str> = results.iter().map(|p| p.transcript_id.as_str()).collect();
        assert_eq!(
            ids,
            vec!["NM_TX1.2"],
            "expected only the curated, current-version transcript, got {:?}",
            ids
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

        /// An `NC_` chromosome parent. Unlike an `NG_`/`LRG_` parent it shares
        /// cdot's coordinate frame, so `project_to_genomic` emits it unchanged
        /// (no re-anchoring, no decline) — the right vehicle for exercising the
        /// c./n./r. → g. coordinate math without tripping the NG_/LRG_ decline
        /// (#655). GRCh38 primary (`.11`) to match the test cdot fixture.
        fn nc_parent() -> Accession {
            Accession::new("NC", "000001", Some(11))
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
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("plus-strand c. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant, got: {}", out),
            };
            // Parent NC_000001.11, plus-strand → g.1003C>A.
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1003C>A"),
                "expected NC_000001.11...:g.1003C>A, got: {}",
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
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("minus-strand c. → g. should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1005G>A"),
                "minus-strand c.4C>T should revcomp to NC_000001.11...:g.1005G>A, got: {}",
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
            let cds = attach_genomic_context_cds(cds, nc_parent());
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

        // -- 3b: minus-strand deletion is re-normalized in genomic space (#737) -

        /// A minus-strand transcript on `NC_000001.11` whose CDS reads
        /// `GTAAAAAAA` (cdna 1..9), so the genome (its reverse complement) carries
        /// a `T` homopolymer. cdot exon `[1000, 1009)` minus-strand maps c.N →
        /// genome `1009 - N`, so the transcript's poly-A run (cdna 3..9) is the
        /// genome poly-T run at g.1000..1006. A single-base deletion anywhere in
        /// that run is the *same* variant; its spec-canonical genomic form is the
        /// 3'-most (highest-coordinate) position, g.1006del.
        ///
        /// Keyed under `NC_000001.11` (not `chr1`) so `normalize` can fetch the
        /// genomic bases back through the same accession the re-anchored output
        /// carries — exercising the #737 genomic-frame renormalization end to end.
        fn make_minus_homopolymer_provider_and_projector() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NM_HOMO_MINUS.1".to_string(),
                CdotTranscript {
                    gene_name: Some("HOMOGENE".to_string()),
                    contig: "NC_000001.11".to_string(),
                    strand: ProvStrand::Minus,
                    exons: vec![[1000, 1009, 0, 9]],
                    cds_start: Some(0),
                    cds_end: Some(9),
                    gene_id: None,
                    protein: Some("NP_HOMO_MINUS.1".to_string()),
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);

            let mut provider = MockProvider::new();
            provider.add_transcript(Transcript {
                id: "NM_HOMO_MINUS.1".to_string(),
                gene_symbol: Some("HOMOGENE".to_string()),
                strand: TxStrand::Minus,
                sequence: Some("GTAAAAAAA".to_string()),
                cds_start: Some(1),
                cds_end: Some(9),
                exons: vec![Exon::new(1, 1, 9)],
                chromosome: Some("NC_000001.11".to_string()),
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
            // Genome forward: prefix + "TTTTTTTAC" (g.1000..1008) + suffix, so the
            // poly-T run is g.1000..1006 (revcomp of the transcript poly-A run).
            let prefix = "N".repeat(999);
            let suffix = "N".repeat(100);
            provider.add_genomic_sequence("NC_000001.11", format!("{}TTTTTTTAC{}", prefix, suffix));
            (projector, provider)
        }

        #[test]
        fn project_to_genomic_minus_strand_deletion_renormalizes_to_genomic_3prime() {
            let (projector, provider) = make_minus_homopolymer_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.9del is the transcript-3'-most representation of the deletion in
            // the poly-A run. On the minus strand that maps to g.1000 — the
            // genome-5' end of the poly-T run. The spec-canonical genomic output
            // is the genome-3'-most position, g.1006del. Before #737 the projector
            // emitted the un-renormalized coordinate image (g.1000del); the
            // re-anchored output is now normalized in its own frame.
            let cds = CdsVariant {
                accession: parse_accession("NM_HOMO_MINUS.1"),
                gene_symbol: Some("HOMOGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(9)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("minus-strand c.9del should project to g.");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant, got: {}", out),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let start = g
                .loc_edit
                .location
                .start
                .inner()
                .expect("start should be concrete");
            assert_eq!(
                start.base, 1006,
                "minus-strand c.9del must re-normalize to the genome-3' end of the \
                 poly-T run (g.1006del), got g.{}del (#737)",
                start.base
            );
        }

        /// Companion to the g.-only test above, exercising the *other* #737
        /// re-anchor site: the multi-axis `.genomic` field built in
        /// `project_single_inner` (via `project_variant`). The two sites apply
        /// the identical genomic-frame renormalization and could drift silently,
        /// so assert the same minus-strand deletion lands at the genome-3' end of
        /// the poly-T run (g.1006del) through a `VariantProjection.genomic`, not
        /// just through the g.-only `project_to_genomic`.
        #[test]
        fn project_variant_minus_strand_deletion_renormalizes_genomic_axis_to_3prime() {
            let (projector, provider) = make_minus_homopolymer_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.9del — the transcript-3'-most representation of the deletion in
            // the poly-A run — maps to the genome-5' end of the poly-T run. The
            // spec-canonical genomic output is the genome-3'-most position,
            // g.1006del; #737 renormalizes the re-anchored `.genomic` axis here.
            let cds = CdsVariant {
                accession: parse_accession("NM_HOMO_MINUS.1"),
                gene_symbol: Some("HOMOGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(9)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_HOMO_MINUS.1")
                .expect("minus-strand c.9del should project (multi-axis)");
            let g = match proj.genomic {
                Some(HgvsVariant::Genome(ref g)) => g,
                other => panic!("expected a Genome `.genomic` axis, got: {:?}", other),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let start = g
                .loc_edit
                .location
                .start
                .inner()
                .expect("start should be concrete");
            assert_eq!(
                start.base, 1006,
                "the multi-axis `.genomic` field must re-normalize to the genome-3' \
                 end of the poly-T run (g.1006del), got g.{}del (#737)",
                start.base
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
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("intronic c.10+5del should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            // Expect a deletion at g.1014 on NC_000001.11.
            assert_eq!(g.accession.to_string(), "NC_000001.11");
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

        /// An `NG_`-genomic input is already in its parent frame, so
        /// `project_to_genomic` must return it unchanged. Before #702 it was fed
        /// through `reanchor_genome_output`, which mistook the parent coordinate
        /// for an `NC_` coordinate and declined (`UnsupportedProjection`: endpoint
        /// outside the placed span) instead of passing the variant through.
        #[test]
        fn project_to_genomic_idempotent_on_ng_genome_input() {
            let (projector, mut provider) = make_nc_two_transcript_setup();
            // NG_900.1 placed on NC_000001.11, plus strand: NG_ base 1 == nc 1000.
            provider.add_genomic_placement(
                "NG_900.1",
                crate::reference::GenomicPlacement {
                    nc: Accession::new("NC", "000001", Some(11)),
                    parent_start: 1,
                    nc_start: 1000,
                    nc_end: 1008,
                    strand: crate::reference::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider);

            // NG_900.1:g.4C>A is in the NG_ frame; base 4 (< nc_start 1000) would
            // be rejected by `nc_to_parent` if treated as an NC_ coordinate.
            let g = GenomeVariant {
                accession: parse_accession("NG_900.1"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Genome(g);
            let out = vp
                .project_to_genomic(&input)
                .expect("NG_ genome input should pass through unchanged");
            assert_eq!(
                out.to_string(),
                input.to_string(),
                "NG_ genome input should be idempotent under project_to_genomic"
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
            let tx = attach_genomic_context_tx(tx, nc_parent());
            let input = HgvsVariant::Tx(tx);
            let out = vp
                .project_to_genomic(&input)
                .expect("n. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1004A>G"),
                "expected NC_000001.11...:g.1004A>G, got: {}",
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
            // DNA-typed bases. RnaPos shape mirrors CdsPos so the coordinate
            // math is identical. (Plus-strand `Base::U` inputs are translated
            // U→T on the g. output — covered separately by
            // `project_to_genomic_rna_u_base_translated_to_t`.)
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // r.4 on NM_TEST.1 (plus strand) → g.1003 on NC_000001.11.
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
            rna.accession = rna.accession.with_genomic_context(nc_parent());
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
                "expected NC_000001.11...:g.1003C>A, got: {}",
                s
            );
            assert_eq!(g.accession.to_string(), "NC_000001.11");
        }

        /// A plus-strand `r.` input carrying `Base::U` is translated to DNA
        /// (`U`→`T`) on the g. output, so the emitted g. variant is valid DNA
        /// rather than an `…U>…` shape (#395 item 4).
        #[test]
        fn project_to_genomic_rna_u_base_translated_to_t() {
            use crate::hgvs::interval::RnaInterval;
            use crate::hgvs::location::RnaPos;
            use crate::hgvs::variant::RnaVariant;

            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // r.4u>g on NM_TEST.1 (plus strand) → g.1003; the RNA `u` reference
            // base must surface as DNA `T`, giving g.1003T>G (ref not validated
            // against the genome, so the edit's own bases are emitted).
            let mut rna = RnaVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    RnaInterval::point(RnaPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::U,
                        alternative: Base::G,
                    },
                ),
            };
            rna.accession = rna.accession.with_genomic_context(nc_parent());
            let out = vp
                .project_to_genomic(&HgvsVariant::Rna(rna))
                .expect("r. → g. projection should succeed");
            let s = out.to_string();
            assert!(
                s.contains(":g.1003T>G"),
                "RNA u> reference should emit DNA T (g.1003T>G), got: {s}"
            );
            assert!(
                !s.contains('U'),
                "g. output must not carry an RNA U base, got: {s}"
            );
        }

        // -- 8b: NG_/LRG_ re-anchor decline → error, not invalid HGVS (#655) ---

        /// An `NG_` parent with no registered placement cannot be re-anchored:
        /// cdot carries only the transcript's chromosome (`NC_`) alignment.
        /// `project_to_genomic` must decline rather than stamp the chromosome
        /// coordinate under the `NG_` accession (invalid HGVS) (#655).
        #[test]
        fn project_to_genomic_ng_parent_without_placement_declines() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
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
            let err = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect_err("NG_ parent without a placement must decline, not emit invalid HGVS");
            match err {
                FerroError::UnsupportedProjection { reason } => assert!(
                    reason.contains("NG_TEST.1") && reason.contains("no chromosomal placement"),
                    "expected a no-placement decline naming NG_TEST.1, got: {reason}"
                ),
                other => panic!("expected UnsupportedProjection, got: {other:?}"),
            }
        }

        /// A placement exists for the `NG_` parent, but the variant's chromosome
        /// coordinate falls outside the placed span, so it cannot be re-anchored
        /// into the parent frame. Decline rather than emit a chromosome
        /// coordinate under the `NG_` accession (#655).
        #[test]
        fn project_to_genomic_ng_parent_endpoint_outside_placement_declines() {
            let (projector, mut provider) = make_test_provider_and_projector();
            // NM_TEST.1 c.4 maps to chromosome 1003; place NG_TEST.1 on a
            // disjoint span [2000, 2010] so nc_to_parent(1003) declines.
            provider.add_genomic_placement(
                "NG_TEST.1",
                crate::reference::GenomicPlacement {
                    nc: Accession::new("NC", "000001", Some(11)),
                    parent_start: 1,
                    nc_start: 2000,
                    nc_end: 2010,
                    strand: crate::reference::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider);
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
            let err = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect_err("an endpoint outside the placed span must decline");
            match err {
                FerroError::UnsupportedProjection { reason } => assert!(
                    reason.contains("outside the placed genomic span"),
                    "expected an out-of-span decline, got: {reason}"
                ),
                other => panic!("expected UnsupportedProjection, got: {other:?}"),
            }
        }

        /// De-anchor path: when the genomic axis cannot be re-anchored into the
        /// `NG_` parent frame (here the endpoint is outside the placement span),
        /// `frame_projection_owned` drops the genomic axis (reports `None`)
        /// rather than stamping the parent accession on an out-of-frame
        /// chromosome coordinate (#655). The transcript-relative coding form is
        /// still re-framed under the parent — degrading the unframable axis is
        /// preferable to dropping the whole projection.
        #[test]
        fn frame_projection_owned_drops_genomic_when_reanchor_declines() {
            let parent = ng_parent("900", 1);
            let placement = crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1010,
                strand: crate::reference::Strand::Plus,
            };
            // Genomic axis at chromosome 5000 — well outside the placed span.
            let genomic = HgvsVariant::Genome(GenomeVariant {
                accession: Accession::new("NC", "000001", Some(11)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(5000)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            });
            let coding = HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_900.1"),
                gene_symbol: Some("GENE900".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            });
            let proj = VariantProjection {
                genomic: Some(genomic),
                coding: Some(coding),
                noncoding: None,
                protein: None,
                rna: None,
                transcript_id: "NM_900.1".to_string(),
                gene_symbol: Some("GENE900".to_string()),
                is_frameshift: false,
                is_intronic: false,
                is_utr: false,
            };
            let framed = VariantProjector::<MockProvider>::frame_projection_owned(
                proj,
                &parent,
                Some(&placement),
            );
            assert!(
                framed.genomic.is_none(),
                "genomic axis should be dropped when it cannot be re-anchored, got: {:?}",
                framed.genomic
            );
            // The coding axis is still re-framed under the NG_ parent.
            let coding_str = framed
                .coding
                .as_ref()
                .expect("coding axis retained")
                .to_string();
            assert!(
                coding_str.starts_with("NG_900.1("),
                "coding should be framed under the NG_ parent, got: {coding_str}"
            );
        }

        /// `reanchor_genome_to_parent` declines an endpoint with no single
        /// resolved position (a compound `(a_b)` `UncertainBoundary::Range`)
        /// with `InvalidCoordinates` rather than stamping the parent accession
        /// on an unresolved coordinate (#655). This decline is defensive — the
        /// public `project_to_genomic` pivot rejects uncertain endpoints
        /// earlier — so it is exercised here directly.
        #[test]
        fn reanchor_genome_to_parent_declines_uncertain_endpoint() {
            use crate::hgvs::interval::UncertainBoundary;
            use crate::hgvs::Mu;
            let placement = crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1010,
                strand: crate::reference::Strand::Plus,
            };
            // The end boundary is a compound range, so `.inner()` is `None`:
            // there is no single position to re-anchor.
            let mut interval = GenomeInterval::new(GenomePos::new(1003), GenomePos::new(1006));
            interval.end = UncertainBoundary::range(
                Mu::Certain(GenomePos::new(1004)),
                Mu::Certain(GenomePos::new(1006)),
            );
            let gv = GenomeVariant {
                accession: ng_parent("900", 1),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    interval,
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let err = VariantProjector::<MockProvider>::reanchor_genome_to_parent(gv, &placement)
                .expect_err("an uncertain endpoint must decline");
            match err {
                FerroError::InvalidCoordinates { msg } => assert!(
                    msg.contains("single resolved position"),
                    "expected an uncertain-position decline, got: {msg}"
                ),
                other => panic!("expected InvalidCoordinates, got: {other:?}"),
            }
        }

        // -- 8c: single-transcript project_variant re-anchors / degrades the
        //        genomic axis for an NG_/LRG_ parent (#702) --------------------

        /// `project_variant` re-anchors the `.genomic` axis into an `NG_`
        /// parent's own frame when a placement is known (#480/#702), rather than
        /// leaving the chromosome-frame pivot stamped with the `NG_` accession.
        #[test]
        fn project_variant_reanchors_genomic_into_ng_parent_with_placement() {
            let (projector, mut provider) = make_test_provider_and_projector();
            // NM_TEST.1 c.4 maps to chromosome 1003; place NG_TEST.1 so chromosome
            // 1000 == NG base 1, hence chromosome 1003 == NG_TEST.1:g.4.
            provider.add_genomic_placement(
                "NG_TEST.1",
                crate::reference::GenomicPlacement {
                    nc: Accession::new("NC", "000001", Some(11)),
                    parent_start: 1,
                    nc_start: 1000,
                    nc_end: 1010,
                    strand: crate::reference::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider);
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
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_TEST.1")
                .expect("project_variant should succeed");
            let g = proj
                .genomic
                .expect("genomic axis should be present and re-anchored");
            let s = g.to_string();
            assert!(
                s.starts_with("NG_TEST.1") && s.contains(":g.4C>A"),
                "expected re-anchored NG_TEST.1:g.4C>A, got: {s}"
            );
            assert!(proj.coding.is_some(), "coding axis should be present");
        }

        /// `project_variant` drops the `.genomic` axis to `None` for an `NG_`
        /// parent with no known placement rather than stamping the `NG_`
        /// accession on a chromosome coordinate (#655/#702). The coding/protein
        /// axes — still valid — are kept; unlike the g.-only `project_to_genomic`,
        /// the multi-axis path degrades rather than hard-errors.
        #[test]
        fn project_variant_drops_genomic_for_ng_parent_without_placement() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
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
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_TEST.1")
                .expect("project_variant should still succeed (multi-axis)");
            assert!(
                proj.genomic.is_none(),
                "genomic axis should be dropped for an unframable NG_ parent, got: {:?}",
                proj.genomic
            );
            assert!(proj.coding.is_some(), "coding axis should be retained");
        }

        /// An `NC_` chromosome parent shares cdot's frame, so `project_variant`
        /// emits its `.genomic` unchanged (regression guard for the #702
        /// re-anchor path).
        #[test]
        fn project_variant_emits_nc_parent_genomic_unchanged() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
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
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_TEST.1")
                .expect("project_variant should succeed");
            let g = proj.genomic.expect("genomic axis present for NC_ parent");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1003C>A"),
                "expected NC_000001.11:g.1003C>A unchanged, got: {s}"
            );
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
            // 999 leading Ns so HGVS g.1000 (the exon's `genome_start`) maps to
            // 0-based index 999, matching every other fixture and the real
            // genome-coordinate convention (HGVS g.X == 0-based X-1). The exon
            // genome and transcript sequences are then identical, so the #644
            // sequence-aware correction is a no-op here.
            let prefix = "N".repeat(999);
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
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);

            let out = vp
                .project_to_genomic(&input)
                .expect("c.*1A>G should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.contains(":g.1012A>G"),
                "expected NC_000001.11...:g.1012A>G for c.*1A>G, got: {}",
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
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);

            let out = vp
                .project_to_genomic(&input)
                .expect("c.-1T>G should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.contains(":g.1002T>G"),
                "expected NC_000001.11...:g.1002T>G for c.-1T>G, got: {}",
                s
            );
        }

        // -- 15: intronic offset on non-coding transcript → unsupported (#332) -

        // -- special positions (pter/qter/cen) in project_to_genomic ----------

        /// #537: a `pter` CDS marker projects to the parent reference's 5'-most
        /// genomic coordinate (`g.1`), not the old special-sentinel rejection.
        /// (PR #534 keeps the *c.-axis* refusal; this covers the g. axis.) The
        /// raw projection is `g.1del`; the normalizer's 3' rule rolls it
        /// downstream.
        #[test]
        fn project_to_genomic_special_pter_projects_to_parent_terminus() {
            let (projector, mut provider) = make_test_provider_and_projector();
            // pter/qter resolve against the parent reference's length.
            provider.add_genomic_sequence("NG_TEST.1", "A".repeat(40));
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
            let g = vp
                .project_to_genomic(&input)
                .expect("pter must project to the parent's 5' terminus");
            assert_eq!(g.to_string(), "NG_TEST.1:g.1del");
        }

        /// A *mixed* numeric/`qter` range (`c.1_qter`) is out of #537 scope —
        /// the numeric endpoint still needs a transcript mapping — so it remains
        /// `UnsupportedProjection`. (A pure `c.qter` resolves to the parent
        /// terminus; that is covered by the `tests/it/issue_537_*` integration
        /// tests, which have a parent length to resolve against.)
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
            // #742: genome bounds load in the HGVS 1-based convention (+1 vs the
            // raw cdot 0-based alt_start/alt_end).
            assert_eq!(tx38.exons, vec![[20001, 20011, 0, 10]]);
            assert_eq!(tx38.cds_start, Some(0));
            assert_eq!(tx38.cds_end, Some(10));
            let tx37 = cdot
                .get_transcript_on_build("NM_MB_TEST.1", "GRCh37")
                .expect("GRCh37 view must be populated in alt_build_transcripts");
            assert_eq!(tx37.contig, "NC_000001.10");
            assert_eq!(tx37.exons, vec![[10001, 10011, 0, 10]]);
        }

        /// c.5A>T against NC_000001.10(NM_MB_TEST.1) → expect GRCh37
        /// coordinates (g.10005; cdot's internal exon[0] is HGVS-numeric,
        /// so c.1 sits at g.10001 and c.5 at g.10005). Pre-fix, this
        /// would have used the primary (GRCh38) view and produced
        /// ~g.20005, despite the parent's GRCh37 version.
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
                        s.contains(":g.10005A>T"),
                        "expected GRCh37 coord g.10005A>T, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// c.5A>T against NC_000001.11(NM_MB_TEST.1) → expect GRCh38
        /// coordinates (g.20005). Sibling to the GRCh37 case; together
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
                        s.contains(":g.20005A>T"),
                        "expected GRCh38 coord g.20005A>T, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// g. → c.: NC_000001.10:g.10005A>T projected onto NM_MB_TEST.1
        /// → expect c.5A>T. Pre-fix, `project_single_inner` used the
        /// primary (GRCh38) `transcript_genome_span` of [20000, 20010),
        /// so position 10005 fell outside and the call returned
        /// `TranscriptNotOverlapping`. Post-fix the GRCh37 span
        /// [10000, 10010) is consulted and the projection succeeds.
        #[test]
        fn project_normalized_uses_grch37_view_for_nc_v10_input() {
            let vp = make_multi_build_projector("GRCh38");
            let g = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10005)),
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

        /// NG_* parents are build-agnostic per the HGVS spec — the chromosome
        /// (`NC_`) pivot falls back to the primary cdot view rather than forcing
        /// an alt-build lookup. With primary=GRCh38 the c.5 position must come
        /// back at the GRCh38 exon's position (g.20005) regardless of any GRCh37
        /// alt-build data.
        ///
        /// This is asserted on the `project_to_genomic_nc` pivot, which is where
        /// the build fallback lives. The public `project_to_genomic` re-anchors
        /// its output and now *declines* an NG_/LRG_ parent with no known
        /// placement rather than stamping a chromosome coordinate under the NG_
        /// accession (#655 — see `project_to_genomic_ng_parent_without_placement_declines`).
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
                .project_to_genomic_nc(&HgvsVariant::Cds(cds))
                .expect("NG-parented variant should pivot via the primary cdot view");
            match out {
                HgvsVariant::Genome(g) => {
                    let s = g.to_string();
                    assert!(
                        s.starts_with("NG_MBTEST.1"),
                        "expected pivot to inherit NG_ parent accession, got: {}",
                        s
                    );
                    assert!(
                        s.contains(":g.20005A>T"),
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

            // c.5del on GRCh37 (g.10005del on NC_000001.10).
            let g37 = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10005)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            vp.project_normalized(&HgvsVariant::Genome(g37), "NM_MB_TEST.1")
                .expect("GRCh37 deletion projection must succeed");

            // c.5del on GRCh38 (g.20005del on NC_000001.11).
            let g38 = GenomeVariant {
                accession: nc_chr1(11),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(20005)),
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
                    GenomeInterval::point(GenomePos::new(20005)),
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
                    GenomeInterval::point(GenomePos::new(10005)),
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
                    GenomeInterval::point(GenomePos::new(20005)),
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

        // ---------------------------------------------------------------------
        // `--assembly` override for build-agnostic NG_ inputs (#715)
        //
        // Builds on the two-build cdot fixture above by giving NM_MB_TEST.1's
        // NG_ parent (NG_000999.1) per-build chromosomal placements with
        // *different* parent offsets, so the re-anchored NG_ coordinate reveals
        // which build's placement was selected: GRCh38 → g.105, GRCh37 → g.205.
        // This is the end-to-end GRCh37 NG_ path #713 deferred (only reachable
        // once the override supplies a build for a bare NG_).
        // ---------------------------------------------------------------------

        const NG_PARENT: &str = "NG_000999.1";

        /// Two-build projector whose NG_000999.1 has a GRCh38 placement
        /// (NC_000001.11, parent_start 100) and — when `with_grch37` — a GRCh37
        /// placement (NC_000001.10, parent_start 200).
        fn make_ng_assembly_projector(with_grch37: bool) -> VariantProjector<MockProvider> {
            use crate::reference::transcript::{
                Exon, GenomeBuild as RefGenomeBuild, ManeStatus, Strand as TxStrand,
                Transcript as RefTranscript,
            };
            use std::sync::OnceLock;

            let cdot =
                CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), "GRCh38")
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
            provider.add_genomic_placement(
                NG_PARENT,
                crate::reference::provider::GenomicPlacement {
                    nc: parse_accession("NC_000001.11"),
                    parent_start: 100,
                    nc_start: 20000,
                    nc_end: 20010,
                    strand: crate::reference::transcript::Strand::Plus,
                },
            );
            if with_grch37 {
                provider.add_genomic_placement(
                    NG_PARENT,
                    crate::reference::provider::GenomicPlacement {
                        nc: parse_accession("NC_000001.10"),
                        parent_start: 200,
                        nc_start: 10000,
                        nc_end: 10010,
                        strand: crate::reference::transcript::Strand::Plus,
                    },
                );
            }
            VariantProjector::new(projector, provider)
        }

        /// c.5A>T on NM_MB_TEST.1 with a *bare* NG_000999.1 parent (no build).
        fn ng_cds_input() -> HgvsVariant {
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_MB_TEST.1")
                    .with_genomic_context(parse_accession(NG_PARENT)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            })
        }

        /// `--assembly GRCh37` makes a bare NG_ select its GRCh37 placement and
        /// re-anchor against the GRCh37 chromosome view (g.205, not g.105).
        #[test]
        fn assembly_override_selects_grch37_placement_for_bare_ng() {
            let vp = make_ng_assembly_projector(true).with_assembly(Some("GRCh37"));
            let out = vp
                .project_to_genomic(&ng_cds_input())
                .expect("GRCh37 override should re-anchor via the GRCh37 placement");
            let s = out.to_string();
            assert!(
                s.contains(":g.205A>T"),
                "expected the GRCh37 placement (parent_start 200 -> g.205), got: {s}"
            );
        }

        /// Without the override a bare NG_ keeps the GRCh38-preferred placement
        /// (g.105) — adding a GRCh37 placement must not shift the default (#713).
        #[test]
        fn no_assembly_override_keeps_grch38_placement_for_bare_ng() {
            let vp = make_ng_assembly_projector(true);
            let out = vp
                .project_to_genomic(&ng_cds_input())
                .expect("default should re-anchor via the GRCh38 placement");
            let s = out.to_string();
            assert!(
                s.contains(":g.105A>T"),
                "expected the GRCh38 placement (parent_start 100 -> g.105), got: {s}"
            );
        }

        /// `--assembly GRCh37` with no GRCh37 placement present declines rather
        /// than mis-anchoring onto the GRCh38 placement (#655/#702 guard).
        #[test]
        fn assembly_override_declines_when_grch37_placement_absent() {
            let vp = make_ng_assembly_projector(false).with_assembly(Some("GRCh37"));
            let err = vp
                .project_to_genomic(&ng_cds_input())
                .expect_err("must decline when the requested build has no placement");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected an UnsupportedProjection decline, got: {err:?}"
            );
        }

        /// A build the input accession already encodes (NC_000001.11) wins over a
        /// contradictory `--assembly GRCh37`: the override only fills in for
        /// build-agnostic inputs, never reinterprets an explicit accession (#715).
        #[test]
        fn input_encoded_build_wins_over_assembly_override() {
            let vp = make_ng_assembly_projector(true).with_assembly(Some("GRCh37"));
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
                .expect("explicit NC_*.11 parent projects on GRCh38 despite the override");
            let s = out.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.20005A>T"),
                "input-encoded GRCh38 (g.20005) must win over --assembly GRCh37 (g.10005), got: {s}"
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
    fn project_nonaug_init_codon_reports_unknown_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_NOATG.1 starts with CTG (non-AUG). c.1C>G affects the initiation
        // codon → unknown consequence, reported as the spec-correct whole-protein
        // p.? (not p.(Met1?), since residue 1 is not Met) without translating
        // (#771). Before the fix this declined via the #625 non-ATG guard.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_NOATG.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("an initiation-codon variant must report p.?");
        assert_eq!(format!("{protein}"), "NP_NOATG.1:p.?");
    }

    #[test]
    fn project_nonaug_downstream_variant_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A NON-initiation variant on the non-AUG transcript must still decline
        // (the #625 guard: translating would risk the wrong frame).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_NOATG.1")
            .expect("projection should succeed");
        assert!(
            proj.protein.is_none(),
            "non-initiation variant on a non-ATG transcript must decline (#625)"
        );
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

    #[test]
    fn project_whole_exon_deletion_predicts_inframe() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.4-1_6+1del: both endpoints intronic, brackets the exonic span c.4_6
        // (codon 2, CGC = Arg). Clamp [4, 6] → in-frame single-residue deletion.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::with_offset(4, -1), CdsPos::with_offset(6, 1)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("whole-exon deletion projection should succeed");
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2del)");
    }

    #[test]
    fn project_pure_intron_deletion_has_no_protein() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.4+1_4+9del: both offsets on base 4 → lo 5 > hi 4 → no exonic CDS → None.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::with_offset(4, 1), CdsPos::with_offset(4, 9)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("pure-intron deletion projection should succeed (no protein)");
        assert!(
            proj.protein.is_none(),
            "pure-intron deletion must not predict a protein"
        );
    }

    #[test]
    fn project_partial_exon_deletion_has_no_protein() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.4-1_5del: one exonic endpoint (c.5) → partial-exon / splice → None.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::with_offset(4, -1), CdsPos::new(5)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("partial-exon deletion projection should succeed (no protein)");
        assert!(
            proj.protein.is_none(),
            "partial-exon deletion must not predict a protein"
        );
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

    // ------------------------------------------------------------------
    // Direct n./r.→CDS→p. path for bare transcript inputs (#506)
    //
    // A bare non-coding (`n.`) or RNA (`r.`) input on a coding transcript
    // (no `genomic_context` parent) carries no genome alignment, so the
    // c.→g.→CDS roundtrip cannot run. But an `n.` base is a 1-based
    // transcript position: convert it to a CDS position via the
    // transcript's `cds_start` (exon/CIGAR-aware `tx_to_cds`), then run the
    // same protein-consequence prediction the coding path uses. The
    // resulting projection has `genomic = None`.
    // ------------------------------------------------------------------

    /// A non-coding transcript (no CDS) for the "no protein on n. input"
    /// guard. Sequence "ACGTACGTAC" (10 nt), exon tx 1..10, plus strand, no
    /// `cds_start`/`cds_end`.
    fn make_noncoding_provider_and_projector() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NR_TEST.1".to_string(),
            CdotTranscript {
                gene_name: Some("NCGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1010, 0, 10]],
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
            gene_symbol: Some("NCGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ACGTACGTAC".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 10)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1009),
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
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ACGTACGTAC", suffix));
        (projector, provider)
    }

    #[test]
    fn project_bare_n_substitution_predicts_protein_without_genome() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1's CDS begins at transcript position 1 (sequence ATGCGCTAA,
        // no 5'UTR), so n.4 == c.4. n.4C>A: codon 2 CGC(Arg) → AGC(Ser).
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare n. input should project via the direct n.→CDS→p. path");
        assert!(
            proj.genomic.is_none(),
            "bare n. input should report genomic = None, got {:?}",
            proj.genomic
        );
        // The non-coding axis is the input itself; the coding axis is derived.
        assert!(
            proj.noncoding.is_some(),
            "noncoding (n.) form should be present"
        );
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2Ser)");
    }

    #[test]
    fn project_bare_n_populates_coding_axis() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // n.4C>A → c.4C>A (cds_start at tx pos 1). The coding form reframes the
        // transcript-relative position to CDS-relative; the edit is unchanged.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare n. input should project");
        let coding = proj
            .coding
            .as_ref()
            .expect("coding (c.) axis should be derived for a coding transcript");
        assert!(
            matches!(coding, HgvsVariant::Cds(_)),
            "coding form should be a c. (Cds) variant, got {coding:?}"
        );
        assert!(
            coding.to_string().contains(":c.4C>A"),
            "n.4C>A should reframe to c.4C>A, got {coding}"
        );
    }

    #[test]
    fn project_bare_r_substitution_predicts_protein_without_genome() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{LocEdit, RnaVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // r.4c>a mirrors n.4C>A on the sense strand. The RNA edit carries
        // `Base::U`-style RNA bases (here `c`/`a`, which are DNA-compatible);
        // the direct path translates U→T before reading the DNA CDS, so the
        // protein consequence matches the coding path: p.(Arg2Ser).
        let variant = HgvsVariant::Rna(RnaVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare r. input should project via the direct r.→CDS→p. path");
        assert!(
            proj.genomic.is_none(),
            "bare r. input should report genomic = None, got {:?}",
            proj.genomic
        );
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2Ser)");
    }

    #[test]
    fn project_bare_r_with_u_base_predicts_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{LocEdit, RnaVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // r.2u>a: tx pos 2 is the `T` of ATG (RNA reads `u`). The U→T
        // translation must run so the ref base matches the DNA CDS. c.2T>A
        // changes the initiation codon ATG→AAG ⇒ p.(Met1?).
        let variant = HgvsVariant::Rna(RnaVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(2)),
                NaEdit::Substitution {
                    reference: Base::U,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare r. input with a U base should project");
        let protein = proj
            .protein
            .expect("a start-codon r. substitution should predict a protein");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Met1?)");
    }

    #[test]
    fn project_bare_n_intronic_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // n.4+5A>G — intronic offset on a bare NM_. Without a genome alignment
        // the intronic position cannot be placed, so protein prediction is not
        // attempted: protein = None, no panic.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::with_offset(4, 5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare n. intronic projection should succeed (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.protein.is_none(),
            "intronic n. variant should not predict a protein consequence"
        );
        assert!(proj.is_intronic, "is_intronic flag should be set");
    }

    #[test]
    fn project_bare_n_on_noncoding_transcript_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_noncoding_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NR_TEST.1 has no CDS, so an n. variant on it has no protein
        // consequence — the projection succeeds with protein = None.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NR_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::T,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NR_TEST.1")
            .expect("bare n. on a non-coding transcript should project (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.protein.is_none(),
            "n. on a non-coding transcript must not predict a protein, got {:?}",
            proj.protein
        );
        // The non-coding axis still carries the input n. form.
        assert!(proj.noncoding.is_some(), "noncoding form should be present");
    }

    #[test]
    fn project_bare_n_utr_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        // NM_UTR5.1: 5 nt 5'UTR then CDS at tx pos 6. n.3 lies in the 5'UTR
        // (c.-3) and does not reach the initiation codon, so no protein.
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_UTR5.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(3)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_UTR5.1")
            .expect("bare n. 5'UTR projection should succeed (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.is_utr,
            "is_utr flag should be set for a 5'UTR n. position"
        );
        assert!(
            proj.protein.is_none(),
            "pure-5'UTR n. position must not predict a protein, got {:?}",
            proj.protein
        );
    }

    #[test]
    fn project_bare_n_transcript_id_mismatch_is_rejected() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // An n. input on NM_TEST.1 projected against a different transcript
        // must be rejected, mirroring the coding path's guard.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_OTHER.1")
            .expect_err("transcript_id mismatch must be rejected");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "transcript_id mismatch should be UnsupportedProjection, got {err:?}"
        );
    }
}
