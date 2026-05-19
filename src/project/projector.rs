//! `VariantProjector` orchestrator.

use crate::data::mapping::MappingInfo;
use crate::data::projection::Projector;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::{CdsInterval, TxInterval};
use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
use crate::hgvs::variant::{
    is_frameshift, AlleleVariant, CdsVariant, HgvsVariant, LocEdit, TxVariant,
};
use crate::normalize::{NormalizeConfig, Normalizer};
use crate::project::accession::parse_accession;
use crate::project::edit::transform_edit_for_strand;
use crate::project::protein::{predict_indel_protein, predict_substitution_protein};
use crate::project::result::VariantProjection;
use crate::reference::ReferenceProvider;

pub struct VariantProjector<P: ReferenceProvider + Clone> {
    projector: Projector,
    provider: P,
    normalizer: Normalizer<P>,
}

impl<P: ReferenceProvider + Clone> VariantProjector<P> {
    pub fn new(projector: Projector, provider: P) -> Self {
        let normalizer = Normalizer::with_config(provider.clone(), NormalizeConfig::default());
        Self {
            projector,
            provider,
            normalizer,
        }
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
        use crate::hgvs::interval::{GenomeInterval, UncertainBoundary};
        use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
        use crate::hgvs::uncertainty::Mu;
        use crate::hgvs::variant::{GenomeVariant, HgvsVariant, LocEdit};
        use crate::reference::Strand as RefStrand;

        // Resolve an `UncertainBoundary<T>` to an owned `T`, distinguishing the
        // two `None` returns of `UncertainBoundary::inner()`:
        //   - `Single(Mu::Unknown)` (parsed `?` sentinel) → `UnsupportedProjection`
        //   - `Range { .. }`        (parsed `(a_b)` range) → `InvalidCoordinates`
        // Without this split, a parsed `c.?` slips through as `InvalidCoordinates`
        // and contradicts the documented contract.
        fn resolve_boundary<T: Copy>(
            boundary: &UncertainBoundary<T>,
            coord_label: &str,
            end_label: &str,
        ) -> Result<T, FerroError> {
            match boundary {
                UncertainBoundary::Single(Mu::Certain(t) | Mu::Uncertain(t)) => Ok(*t),
                UncertainBoundary::Single(Mu::Unknown) => Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "project_to_genomic cannot resolve unknown position (`?`) at \
                         {coord_label} {end_label}"
                    ),
                }),
                UncertainBoundary::Range { .. } => Err(FerroError::InvalidCoordinates {
                    msg: format!("{coord_label} interval {end_label} is a range"),
                }),
            }
        }

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
                let s = resolve_boundary(&v.loc_edit.location.start, "c.", "start")?;
                let e = resolve_boundary(&v.loc_edit.location.end, "c.", "end")?;
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    s,
                    e,
                )
            }
            HgvsVariant::Tx(v) => {
                let s = resolve_boundary(&v.loc_edit.location.start, "n.", "start")?;
                let e = resolve_boundary(&v.loc_edit.location.end, "n.", "end")?;
                let to_cds = |p: TxPos| CdsPos {
                    base: p.base,
                    offset: p.offset,
                    utr3: p.downstream,
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
                let s = resolve_boundary(&v.loc_edit.location.start, "r.", "start")?;
                let e = resolve_boundary(&v.loc_edit.location.end, "r.", "end")?;
                // RnaPos shape mirrors CdsPos exactly.
                let to_cds = |p: crate::hgvs::location::RnaPos| CdsPos {
                    base: p.base,
                    offset: p.offset,
                    utr3: p.utr3,
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
        if start_cds.is_unknown() || end_cds.is_unknown() {
            return Err(FerroError::UnsupportedProjection {
                reason: "project_to_genomic cannot resolve `?` position sentinels".to_string(),
            });
        }

        // 3. Require an explicit parent NG/NC accession via `genomic_context`.
        //    Per design (#327) we do NOT synthesize a parent from cdot.
        let parent = accession
            .genomic_context
            .as_deref()
            .cloned()
            .ok_or_else(|| FerroError::UnsupportedProjection {
                reason: format!(
                    "input variant has no parent reference (genomic_context) on accession {}; \
                     project_to_genomic requires an explicit NG/NC parent (see #327)",
                    accession.full()
                ),
            })?;

        // 4. Resolve the transcript via the cdot mapper (the same backing
        //    store used by the g. → c./n. path).  Working directly off cdot
        //    keeps us symmetric with `project_single_inner` and avoids a
        //    second `provider.get_transcript` load whose exon entries may
        //    lack genomic coordinates.
        let transcript_id = accession.transcript_accession();
        let cdot_mapper = self.projector.mapper();
        let cdot_tx = cdot_mapper
            .cdot()
            .get_transcript(&transcript_id)
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
                let result = cdot_mapper.cds_to_genome(&transcript_id, &p)?;
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
        let (contig, pos) = extract_contig_and_pos(variant)?;

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
            let gene_symbol = self
                .projector
                .mapper()
                .cdot()
                .get_transcript(transcript_id)
                .ok_or_else(|| FerroError::ReferenceNotFound {
                    id: transcript_id.to_string(),
                })?
                .gene_name
                .clone();
            return Ok(VariantProjection {
                genomic: original.clone(),
                coding: None,
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

        Ok(VariantProjection {
            genomic: original.clone(),
            coding,
            protein,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift,
            is_intronic,
            is_utr,
        })
    }

    /// Project a single (non-allele) g. variant, assuming it has already been
    /// normalized.
    fn project_single_inner(
        &self,
        normalized: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Require a g. variant.
        let genome_variant = match normalized {
            HgvsVariant::Genome(g) => g.clone(),
            _ => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "VariantProjector currently only accepts g. variants".to_string(),
                });
            }
        };

        let edit = genome_variant
            .loc_edit
            .edit
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::UnsupportedProjection {
                reason: "g. variant has no concrete edit".to_string(),
            })?;

        // Look up the transcript in the cdot mapper.
        let cdot_tx = self
            .projector
            .mapper()
            .cdot()
            .get_transcript(transcript_id)
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
        let normalized_str = normalized.to_string();

        // Compute the overall genomic extent of the transcript from its exons.
        // Exon format: [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end].
        let (tx_genome_start, tx_genome_end) = {
            let exons = &cdot_tx.exons;
            if exons.is_empty() {
                return Err(FerroError::ReferenceNotFound {
                    id: transcript_id.to_string(),
                });
            }
            let starts = exons.iter().map(|e| e[0]).min().unwrap();
            let ends = exons.iter().map(|e| e[1]).max().unwrap();
            (starts, ends)
        };

        // Helper: map one GenomePos → CdsPos, converting out-of-range errors.
        let map_position = |gp: &GenomePos| -> Result<(CdsPos, MappingInfo), FerroError> {
            if gp.base < tx_genome_start || gp.base >= tx_genome_end {
                return Err(FerroError::TranscriptNotOverlapping {
                    variant: normalized_str.clone(),
                    transcript_id: transcript_id.to_string(),
                });
            }
            match mapper.genome_to_cds(transcript_id, gp) {
                Ok(res) => Ok((res.variant, res.info)),
                Err(FerroError::InvalidCoordinates { .. })
                | Err(FerroError::ConversionError { .. }) => {
                    Err(FerroError::TranscriptNotOverlapping {
                        variant: normalized_str.clone(),
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

        // Build the c./n. HGVS variant.
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

        // Derive per-position flags from MappingInfo.
        let is_intronic = info_start.is_intronic || info_end.is_intronic;
        let is_utr = !is_intronic
            && (info_start.in_5utr || info_start.in_3utr || info_end.in_5utr || info_end.in_3utr);

        // Predict protein consequence for CDS variants.
        let mut protein = None;
        if !is_intronic && !is_utr && is_coding {
            // Resolve the protein accession in priority order:
            //   1. `cdot_protein` — explicit value from cdot.protein or from the
            //      ferro reference JSON's `protein_id` field (which the cdot
            //      mapper builders plumb into `CdotTranscript.protein`).
            //   2. RefSeq inference: `NM_*` → `NP_*`, `XM_*` → `XP_*`. Substring
            //      stripping (not replace) avoids mangling accessions whose
            //      suffix happens to contain `NM_`.
            //   3. The transcript ID itself, used as the `p.` accession prefix.
            //      Keeps protein prediction live for references whose transcript
            //      IDs do not match a RefSeq convention (e.g. `MY_GENE-gene.1`
            //      from `ferro convert-gff`). The HGVS protein grammar requires
            //      a `sequence_identifier`, so accession-less `p.` is not a valid
            //      alternative. See #310.
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
            if let Some(prot_acc) = prot_acc {
                match &c_edit {
                    NaEdit::Substitution { .. } => {
                        let tx_for_codon = self.provider.get_transcript(transcript_id)?;
                        protein = Some(predict_substitution_protein(
                            &tx_for_codon,
                            cds_start.base,
                            &c_edit,
                            &prot_acc,
                        )?);
                    }
                    NaEdit::Deletion { .. }
                    | NaEdit::Insertion { .. }
                    | NaEdit::Duplication { .. }
                    | NaEdit::Delins { .. }
                    | NaEdit::Inversion { .. }
                        // Only predict when both CDS positions are concrete exonic positions
                        // (no intronic offsets — already guarded above).
                        if cds_start.offset.is_none()
                            && cds_end.offset.is_none()
                            && cds_start.base > 0
                            && cds_end.base > 0 =>
                    {
                        let tx_for_codon = self.provider.get_transcript(transcript_id)?;
                        match predict_indel_protein(
                            &tx_for_codon,
                            cds_start.base,
                            cds_end.base,
                            &c_edit,
                            &prot_acc,
                        ) {
                            Ok(pv) => protein = Some(pv),
                            // Non-fatal: unsupported edits or missing sequence → leave protein=None.
                            Err(FerroError::UnsupportedProjection { .. })
                            | Err(FerroError::ProteinSequenceUnavailable { .. }) => {}
                            Err(other) => return Err(other),
                        }
                    }
                    _ => {}
                }
            }
        }

        let frameshift = is_frameshift(&coding);
        Ok(VariantProjection {
            genomic: normalized.clone(),
            coding: Some(coding),
            protein,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
        })
    }
}

/// Extract the contig name and a representative 1-based genomic position from
/// an already-normalized `HgvsVariant`.
///
/// For `HgvsVariant::Allele`, the first inner g. variant is used.
///
/// # Contig name resolution
///
/// The contig name is taken directly from the variant's `Accession`. For
/// standard RefSeq genomic accessions (e.g. `NC_000001.11`) the cdot mapper
/// stores contigs under those same names and also aliases UCSC names
/// (`chr1`) to them via `populate_contig_aliases`.  For assembly-notation
/// accessions (`GRCh38(chr1)`) the `chromosome` field of `Accession` is
/// used instead.  Callers that use non-standard accession formats (e.g.
/// plain `chr1` keys) should ensure their cdot data was loaded with matching
/// contig keys.
fn extract_contig_and_pos(variant: &HgvsVariant) -> Result<(String, u64), FerroError> {
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
            // Prefer the chromosome field for assembly-notation accessions.
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
        _ => Err(FerroError::UnsupportedProjection {
            reason: "project_all currently only accepts g. variants".to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::{CdotMapper, CdotTranscript};
    use crate::data::projection::Projector;
    use crate::hgvs::variant::{AllelePhase, AlleleVariant};
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
}
