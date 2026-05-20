//! Normalization engine
//!
//! Implements the core HGVS variant normalization algorithm including
//! 3'/5' shifting and boundary detection.
//!
//! # Coordinate Systems
//!
//! This module uses multiple coordinate systems:
//!
//! | Context | Basis | Type/Notes |
//! |---------|-------|------------|
//! | HGVS variant positions | 1-based | `u64` for genomic, `i64` for CDS/Tx |
//! | Array indexing | 0-based | `usize` for sequence slicing |
//! | Boundaries struct | 1-based | `(start, end)` inclusive |
//! | Shuffle input/output | 0-based | Uses array indices |
//! | Relative positions | 1-based | Positions within fetched window |
//!
//! Key conversions:
//! - `hgvs_pos_to_index(pos)` converts 1-based HGVS position to 0-based index
//! - `index_to_hgvs_pos(idx)` converts 0-based index to 1-based HGVS position
//! - `pos.saturating_sub(1)` manually converts 1-based to 0-based
//! - `idx + 1` manually converts 0-based to 1-based
//!
//! For type-safe coordinate handling, see [`crate::coords`].

pub mod boundary;
pub mod config;
pub(crate) mod merge;
mod overlap;
pub mod rules;
pub mod shuffle;
pub mod validate;

use crate::coords::{hgvs_pos_to_index, index_to_hgvs_pos};
use crate::error::FerroError;
use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::Interval;
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    AllelePhase, AlleleVariant, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant,
    RnaVariant, TxVariant,
};
use crate::hgvs::HgvsVariant as HV;
use crate::reference::transcript::Strand;
use crate::reference::ReferenceProvider;
use boundary::Boundaries;
pub use config::{NormalizeConfig, ShuffleDirection};
use rules::{
    canonicalize_conversion_to_delins, canonicalize_edit, needs_normalization, should_canonicalize,
    DelinsSubedit,
};
use shuffle::shuffle;

/// Check if a CDS position has an unknown (?) offset sentinel value
fn has_unknown_offset_cds(pos: &CdsPos) -> bool {
    matches!(
        pos.offset,
        Some(OFFSET_UNKNOWN_POSITIVE) | Some(OFFSET_UNKNOWN_NEGATIVE)
    )
}

/// Check if a TxPos has an unknown (?) offset sentinel value
fn has_unknown_offset_tx(pos: &TxPos) -> bool {
    matches!(
        pos.offset,
        Some(OFFSET_UNKNOWN_POSITIVE) | Some(OFFSET_UNKNOWN_NEGATIVE)
    )
}

/// Warning generated during normalization.
///
/// This enum is open-ended: each variant owns the fields its code needs.
/// Future warning codes add new variants without touching existing emit sites.
#[derive(Debug, Clone)]
pub enum NormalizationWarning {
    /// Reference sequence mismatch. Stated ref bases in the HGVS expression
    /// do not match the actual reference sequence. Code: `REFSEQ_MISMATCH`.
    RefSeqMismatch {
        /// Human-readable description
        message: String,
        /// What the input claimed as reference
        stated_ref: String,
        /// What the actual reference sequence has
        actual_ref: String,
        /// Position info
        position: String,
        /// Whether the mismatch was actually corrected by the normalizer.
        ///
        /// `true` when the canonical Display drops or rewrites the stated
        /// bases (sub / del / dup / inv: the explicit `sequence` field is
        /// stripped during canonicalization, so the wrong stated bases
        /// never reach the output).
        ///
        /// `false` when the description passes through verbatim — the
        /// normalizer surfaces the warning but cannot rewrite the user's
        /// declared form. This covers `Repeat` and `MultiRepeat`
        /// consistency mismatches (issues #214 / #279): the per-unit
        /// declaration is part of the user's variant description, and
        /// the normalizer declines to second-guess it. (Issue #280.)
        corrected: bool,
    },

    /// Two or more cis-allele edits share identical reference bounds.
    /// The HGVS spec does not define a canonical form for this case;
    /// ferro preserves the input verbatim and emits this warning.
    /// Code: `OVERLAP_CONFLICTING_EDITS`.
    OverlapConflict {
        /// Human-readable description
        message: String,
        /// Accession of the reference sequence
        accession: String,
        /// Coordinate system: "g" | "c" | "n" | "r" | "m"
        coordinate_system: String,
        /// Canonical span text, e.g. "100" or "100_103"
        location: String,
        /// Edit kinds, e.g. ["sub", "sub"] or ["del", "inv"]
        edit_kinds: Vec<String>,
    },

    /// `apply_canonical_split` was unable to canonicalize because the
    /// reference window returned by the provider did not span the same
    /// number of bytes as the HGVS interval. Per HGVS spec
    /// `recommendations/background/refseq.md` §43, this means the
    /// variant is not entirely encompassed by the reference sequence —
    /// strict mode promotes this warning to
    /// `FerroError::VariantExceedsReference`. The variant is returned
    /// unchanged in lenient/silent modes.
    /// Closes-after: #354, #355. Code: `CANONICAL_SPLIT_SKIPPED`.
    CanonicalSplitSkipped {
        /// Human-readable description
        message: String,
        /// Accession of the reference sequence
        accession: String,
        /// HGVS span start (1-based inclusive). Carried so strict-mode
        /// promotion can build a `FerroError::VariantExceedsReference`
        /// with the same span information.
        hgvs_start: u64,
        /// HGVS span end (1-based inclusive).
        hgvs_end: u64,
        /// Number of bytes the HGVS span demanded.
        expected_span: usize,
        /// Number of bytes the provider returned.
        actual_bytes: usize,
    },

    /// A `c.` variant whose start and end positions sit in different
    /// coordinate sub-axes (5'UTR / CDS / 3'UTR). The 3'-rule shuffle has
    /// no well-defined semantics across an axis boundary, so ferro
    /// preserves the canonical input position and emits this warning.
    /// Closes-after: #350. Code: `CROSS_AXIS_VARIANT_NOT_SHUFFLED`.
    CrossAxisVariantNotShuffled {
        /// Human-readable description
        message: String,
        /// Accession of the reference sequence
        accession: String,
        /// Axis of the start position: "5utr" | "cds" | "3utr"
        start_axis: String,
        /// Axis of the end position: "5utr" | "cds" | "3utr"
        end_axis: String,
    },

    /// A 3'-rule shuffle would have crossed a CDS↔UTR coordinate sub-axis
    /// boundary, but the axis clamp constrained the result to the
    /// boundary. Closes-after: #349. Code: `AXIS_CLAMP_APPLIED`.
    AxisClampApplied {
        /// Human-readable description
        message: String,
        /// Accession of the reference sequence
        accession: String,
        /// Shuffle direction that was clamped: "5prime" | "3prime"
        direction: String,
        /// Axis bound that did the clamping: "5utr" | "cds_start" | "cds_end" | "3utr"
        clamp_kind: String,
    },
}

impl NormalizationWarning {
    /// The warning's user-facing code string.
    pub fn code(&self) -> &'static str {
        match self {
            Self::RefSeqMismatch { .. } => "REFSEQ_MISMATCH",
            Self::OverlapConflict { .. } => "OVERLAP_CONFLICTING_EDITS",
            Self::CanonicalSplitSkipped { .. } => "CANONICAL_SPLIT_SKIPPED",
            Self::CrossAxisVariantNotShuffled { .. } => "CROSS_AXIS_VARIANT_NOT_SHUFFLED",
            Self::AxisClampApplied { .. } => "AXIS_CLAMP_APPLIED",
        }
    }

    /// Human-readable message for the warning.
    pub fn message(&self) -> &str {
        match self {
            Self::RefSeqMismatch { message, .. } => message,
            Self::OverlapConflict { message, .. } => message,
            Self::CanonicalSplitSkipped { message, .. } => message,
            Self::CrossAxisVariantNotShuffled { message, .. } => message,
            Self::AxisClampApplied { message, .. } => message,
        }
    }
}

/// Result of normalization with optional warnings
#[derive(Debug, Clone)]
pub struct NormalizeResultWithWarnings {
    /// The normalized variant
    pub result: HgvsVariant,
    /// Warnings generated during normalization
    pub warnings: Vec<NormalizationWarning>,
}

impl NormalizeResultWithWarnings {
    /// Create a new result without warnings
    pub fn new(result: HgvsVariant) -> Self {
        Self {
            result,
            warnings: vec![],
        }
    }

    /// Create a result with warnings
    pub fn with_warnings(result: HgvsVariant, warnings: Vec<NormalizationWarning>) -> Self {
        Self { result, warnings }
    }

    /// Add a warning to the result
    pub fn add_warning(&mut self, warning: NormalizationWarning) {
        self.warnings.push(warning);
    }

    /// Check if there are any warnings
    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }

    /// Check if there's a reference mismatch warning
    pub fn has_ref_mismatch(&self) -> bool {
        self.warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }))
    }
}

/// Main normalizer struct
pub struct Normalizer<P: ReferenceProvider> {
    provider: P,
    config: NormalizeConfig,
}

impl<P: ReferenceProvider> Normalizer<P> {
    /// Create a new normalizer with the given reference provider
    pub fn new(provider: P) -> Self {
        Self {
            provider,
            config: NormalizeConfig::default(),
        }
    }

    /// Create a normalizer with custom configuration
    pub fn with_config(provider: P, config: NormalizeConfig) -> Self {
        Self { provider, config }
    }

    /// Get the configuration
    pub fn config(&self) -> &NormalizeConfig {
        &self.config
    }

    /// Normalize a variant
    ///
    /// In strict mode (default), rejects variants with reference mismatches.
    /// Use `normalize_with_warnings` for lenient mode that corrects mismatches.
    pub fn normalize(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        let result = self.normalize_with_warnings(variant)?;

        // In strict mode, reject if there were reference mismatches.
        if self.config.should_reject_ref_mismatch() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::RefSeqMismatch {
                    position,
                    stated_ref,
                    actual_ref,
                    ..
                } => Some(FerroError::ReferenceMismatch {
                    location: position.clone(),
                    expected: stated_ref.clone(),
                    found: actual_ref.clone(),
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        // Strict mode also rejects W5003 VariantExceedsReference per
        // HGVS spec refseq.md §43 — the variant must be entirely
        // encompassed by the selected reference. Promotes the
        // `CanonicalSplitSkipped` warning to a typed error. Closes
        // #355; matches biocommons hgvs which raises
        // `HGVSInvalidVariantError` for this shape.
        if self.config.should_reject_variant_exceeds_reference() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::CanonicalSplitSkipped {
                    accession,
                    hgvs_start,
                    hgvs_end,
                    expected_span,
                    actual_bytes,
                    ..
                } => Some(FerroError::VariantExceedsReference {
                    accession: accession.clone(),
                    hgvs_start: *hgvs_start,
                    hgvs_end: *hgvs_end,
                    expected_span: *expected_span as u64,
                    actual_bytes: *actual_bytes as u64,
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        Ok(result.result)
    }

    /// Normalize a variant with detailed warnings
    ///
    /// Returns the normalized variant along with any warnings generated during
    /// normalization (e.g., reference sequence mismatches that were auto-corrected).
    /// Use this method when you want to track what corrections were made.
    pub fn normalize_with_warnings(
        &self,
        variant: &HgvsVariant,
    ) -> Result<NormalizeResultWithWarnings, FerroError> {
        let (result, warnings) = match variant {
            HV::Genome(v) => self.normalize_genome(v)?,
            HV::Cds(v) => self.normalize_cds(v)?,
            HV::Tx(v) => self.normalize_tx(v)?,
            HV::Protein(v) => self.normalize_protein(v)?,
            HV::Rna(v) => self.normalize_rna(v)?,
            HV::Mt(v) => self.normalize_mt(v)?,
            HV::Allele(a) => self.normalize_allele(a)?,
            // Circular variants normalize like genomic variants
            HV::Circular(v) => (
                HV::Circular(crate::hgvs::variant::CircularVariant {
                    accession: v.accession.clone(),
                    gene_symbol: v.gene_symbol.clone(),
                    loc_edit: v.loc_edit.clone(),
                }),
                vec![],
            ),
            // RNA fusions pass through unchanged (no normalization needed for fusions)
            HV::RnaFusion(v) => (HV::RnaFusion(v.clone()), vec![]),
            // Null and unknown allele markers pass through unchanged
            HV::NullAllele => (HV::NullAllele, vec![]),
            HV::UnknownAllele => (HV::UnknownAllele, vec![]),
        };

        Ok(NormalizeResultWithWarnings::with_warnings(result, warnings))
    }

    /// Normalize an allele (compound) variant
    ///
    /// Normalizes each variant in the allele individually, with overlap prevention.
    /// After normalization, checks if variants would overlap and constrains shifting
    /// to prevent collisions.
    fn normalize_allele(
        &self,
        allele: &crate::hgvs::variant::AlleleVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Merge first, then normalize each result through the full per-variant
        // pipeline. Pre-normalizing each bracket entry would shift it
        // independently of its siblings, which can collapse adjacent edits
        // onto the same 3'-end position (e.g. `[260delA;261delA]` →
        // `[264del;264del]`) and defeat the strict-adjacency merge. Issue #180.
        //
        // Order:
        //   1. merge raw bracket entries by positional adjacency
        //   2. cis-only: decompose merged delins into [..., inv, ...] (#160)
        //   3. run the full per-variant pipeline on every result — this is
        //      what applies the HGVS 3' rule to the merged anchor (#161, #180)
        //   4. detect post-shift overlaps and emit warnings
        let mut all_warnings = Vec::new();
        let original_len = allele.variants.len();
        let merged_raw =
            merge::merge_consecutive_edits(allele.variants.clone(), allele.phase, &self.provider);

        // Issue #160 + #165: any merged delins (or pre-existing delins
        // that survived merge unchanged) may decompose into a sequence of
        // higher-priority forms per `general.md:56` — `[..., inv, ...]`
        // when an inv-eligible sub-span is present (#160) and/or into
        // separate substitutions when interior positions match the
        // reference (#165). Run the split per merged variant; the helper
        // is a no-op for non-Delins variants and for Delins with nothing
        // to split out. Only applies in cis phase — trans alleles aren't
        // collapsible in the first place.
        let merged_split: Vec<HgvsVariant> =
            if allele.phase == crate::hgvs::variant::AllelePhase::Cis {
                merged_raw
                    .into_iter()
                    .flat_map(|v| self.canonical_split_for_variant(v))
                    .collect()
            } else {
                merged_raw
            };

        // Per-variant pipeline on every merged result. This is the single
        // canonical place where the 3' rule, ins→dup canonicalization, ref
        // validation, etc. apply — a merged variant is semantically a new
        // variant and goes through the same pipeline as any direct input.
        let mut normalized: Vec<HgvsVariant> = Vec::with_capacity(merged_split.len());
        for v in merged_split {
            let r = self.normalize_with_warnings(&v)?;
            all_warnings.extend(r.warnings);
            normalized.push(r.result);
        }

        // Overlap detection runs post-shift so collisions caused by the
        // 3' shift surface alongside input-time ones. Overlap *prevention*
        // is structural — the merge-first ordering above plus the strict
        // `prev.end + 1 == next.start` adjacency check in
        // `merge_consecutive_edits` make it impossible for the normalizer
        // to emit overlapping ranges from non-overlapping inputs.
        all_warnings.extend(crate::normalize::overlap::detect_overlap_conflicts(
            &normalized,
            allele.phase,
        ));

        // HGVS requires consecutive edits in cis to render as a single edit.
        // Only unwrap when a merge actually collapsed multiple sub-variants —
        // pre-existing singleton alleles must round-trip with the Allele
        // wrapper intact for programmatic callers (Display already renders
        // singletons in bare form regardless).
        let result = if allele.phase == crate::hgvs::variant::AllelePhase::Cis
            && original_len > 1
            && normalized.len() == 1
        {
            normalized.pop().unwrap()
        } else {
            HgvsVariant::Allele(crate::hgvs::variant::AlleleVariant::new(
                normalized,
                allele.phase,
            ))
        };

        Ok((result, all_warnings))
    }

    /// Normalize a genomic variant
    fn normalize_genome(
        &self,
        variant: &GenomeVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Genome(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins` before any further work.
        // Pure-syntax; no reference data needed.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = GenomeVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return Ok((HV::Genome(new_variant), vec![]));
        }

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Genome(variant.clone()), vec![]));
        }

        // Get sequence from provider - can't normalize with unknown positions
        let accession = variant.accession.transcript_accession();
        let start = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos.base,
            None => {
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    vec![],
                ))
            }
        };
        let end = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos.base,
            None => {
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    vec![],
                ))
            }
        };

        // Try to get transcript/sequence, fall back to minimal notation if not found
        let window_start = start.saturating_sub(self.config.window_size);
        let seq_result = self.provider.get_sequence(
            &accession,
            window_start,
            end.saturating_add(self.config.window_size),
        );

        let ref_seq = match seq_result {
            Ok(s) => s,
            // Can't do full normalization without sequence, but apply minimal notation
            Err(_) => {
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    vec![],
                ))
            }
        };

        // Adjust coordinates to be relative to the window
        let rel_start = start - window_start;
        let rel_end = end - window_start;

        // Perform normalization
        let (new_rel_start, new_rel_end, new_edit, mut warnings) = self.normalize_na_edit(
            ref_seq.as_bytes(),
            edit,
            rel_start,
            rel_end,
            &Boundaries::new(0, ref_seq.len() as u64),
            false, // genomic context: codon-frame gate does not apply
        )?;

        // Adjust back to genomic coordinates
        let new_start = new_rel_start + window_start;
        let new_end = new_rel_end + window_start;

        // Reconstruct variant with new position (using adjusted coordinates)
        let new_variant = GenomeVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(new_start), GenomePos::new(new_end)),
                new_edit,
            ),
        };

        // Issue #160 + #165: a normalized Delins may decompose under the
        // spec's edit-priority rule (`general.md:56`) — into `[..., inv,
        // ...]` for rev-comp sub-spans and/or into separate subs across
        // interior identities. Returns the variant unchanged for
        // non-Delins or no-decomposition cases.
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Genome(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Normalize a CDS variant
    fn normalize_cds(
        &self,
        variant: &CdsVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Cds(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins`.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = CdsVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return Ok((HV::Cds(new_variant), vec![]));
        }

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Cds(variant.clone()), vec![]));
        }

        // Check for intronic variants or unknown positions - can't normalize those
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Can't normalize variants with unknown (?) offsets - return unchanged
        if has_unknown_offset_cds(start_pos) || has_unknown_offset_cds(end_pos) {
            return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![]));
        }

        // Try to get transcript first - we need it for intronic normalization too.
        //
        // For intronic / boundary-spanning positions we route through
        // `get_transcript_for_variant` so providers that hold cdot data
        // (e.g. `MultiFastaProvider`) can pick a genome build informed by
        // the input's `genomic_context` (NG_* / NC_* parent). For the
        // pure-exonic CDS path the build hint is irrelevant — `get_transcript`
        // is sufficient.
        let accession = variant.accession.transcript_accession();
        let transcript_for_intronic =
            || -> Result<crate::reference::transcript::Transcript, FerroError> {
                self.provider
                    .get_transcript_for_variant(&HV::Cds(variant.clone()))
            };
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            // Can't do full normalization without transcript, but apply minimal notation
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Handle intronic variants specially
        if start_pos.is_intronic() || end_pos.is_intronic() {
            // Switch to the variant-aware lookup so an NG/NC-parented input
            // gets the build-correct chromosome. If the variant-aware lookup
            // fails, fall back to the plain transcript we already fetched.
            let transcript = transcript_for_intronic().unwrap_or(transcript);
            // Check if both positions are intronic and in the same intron
            if start_pos.is_intronic() && end_pos.is_intronic() {
                return self.normalize_intronic_cds(variant, &transcript, start_pos, end_pos, edit);
            }
            // Variant spans exon-intron boundary - normalize in genomic space
            return self.normalize_boundary_spanning_cds(
                variant,
                &transcript,
                start_pos,
                end_pos,
                edit,
            );
        }

        // Convert CDS to transcript coordinates for normalization
        let cds_start = match transcript.cds_start {
            Some(s) => s,
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Calculate transcript positions - return unchanged if position is out of range
        let tx_start = match self.cds_to_tx_pos(start_pos, cds_start, transcript.cds_end) {
            Ok(v) => v,
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        let tx_end = match self.cds_to_tx_pos(end_pos, cds_start, transcript.cds_end) {
            Ok(v) => v,
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Get boundaries. Both cross_boundaries modes route through
        // `get_cds_boundaries_with_axis_info` so the CDS↔UTR axis bound
        // applies regardless of cross mode (closes #337). The
        // exon-vs-full-tx dimension still toggles on
        // `config.cross_boundaries`. We additionally use the un-clamped
        // exon bound to detect:
        //
        //   - Cross-axis variants (#350): `tx_start` and `tx_end` map to
        //     different axes (5'UTR / CDS / 3'UTR). The 3'-rule shuffle
        //     has no well-defined semantics across an axis boundary, so
        //     we preserve the canonical input position and emit
        //     `CrossAxisVariantNotShuffled`.
        //
        //   - Axis-clamp activations (#349): after shuffling, the result
        //     position rests at the axis boundary AND the axis bound is
        //     tighter than the exon bound on that side. Emit
        //     `AxisClampApplied` so callers can flag for human review.
        let axis_info = match boundary::get_cds_boundaries_with_axis_info(
            &transcript,
            tx_start,
            &self.config,
        ) {
            Ok(b) => b,
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        let boundaries = axis_info.clamped.clone();
        let exon_only = axis_info.exon.clone();
        let start_axis = axis_info.axis_region;
        let end_axis = boundary::axis_region_of(&transcript, tx_end);

        // #350: bail on cross-axis variants. Both endpoints must live in
        // the same axis region for a 3'-rule shuffle to be defined.
        if start_axis != end_axis
            && !matches!(start_axis, boundary::AxisRegion::None)
            && !matches!(end_axis, boundary::AxisRegion::None)
        {
            let acc = variant.accession.transcript_accession();
            let warning = NormalizationWarning::CrossAxisVariantNotShuffled {
                message: format!(
                    "{}:c.{}: range straddles {} and {} axes; \
                     3'-rule shuffle is undefined across a CDS\u{2194}UTR \
                     boundary, returning input unchanged",
                    acc,
                    variant.loc_edit.location,
                    start_axis.label(),
                    end_axis.label(),
                ),
                accession: acc,
                start_axis: start_axis.label().to_string(),
                end_axis: end_axis.label().to_string(),
            };
            return Ok((
                HV::Cds(self.canonicalize_cds_variant(variant)),
                vec![warning],
            ));
        }

        // Perform normalization on transcript sequence (CDS context).
        // Coordinate-only transcripts (no cached bases) fall back to the
        // canonicalize-only path, matching the other early-return branches.
        let seq = match transcript.sequence.as_deref() {
            Some(s) => s.as_bytes(),
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        // HGVS spec (repeated.md): the codon-frame restriction
        // (`unit_len % 3 == 0` for repeat notation in `c.` context)
        // applies only to bases inside the CDS proper. UTR positions
        // are exempt:
        //   > This restriction only applies to the coding sequence,
        //   > which does not include the introns or the UTR sequence.
        // The variant is entirely within the CDS iff its transcript-
        // frame span sits between `cds_start` and `cds_end` (inclusive).
        // 5' UTR (`c.-N`) maps to `tx_start < cds_start`; 3' UTR
        // (`c.*N`) maps to `tx_start > cds_end`. Pass `is_coding=false`
        // for any variant whose footprint touches UTR. If `cds_end` is
        // unset we cannot verify the variant lies within CDS proper, so
        // we conservatively treat it as UTR-touching and skip the gate.
        let is_coding = match transcript.cds_end {
            Some(cds_end) => tx_start >= cds_start && tx_end <= cds_end,
            None => false,
        };
        let (new_tx_start, new_tx_end, new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, is_coding)?;

        // #349: detect whether the axis clamp was operative for this
        // shuffle. For 5'-direction shuffles the clamp fires when the
        // result is anchored at `boundaries.left` AND the axis bound is
        // strictly tighter than the exon bound on the left. Symmetric
        // logic for 3'-direction. Skip when there's no axis sub-region
        // (non-coding transcripts).
        //
        // The cheap landed-at-boundary check alone is not sufficient — it
        // can fire even when the unclamped shuffle would have stopped at
        // the same axis-boundary position anyway (e.g. when the reference
        // base immediately past the boundary does not match, so the
        // shuffle was never going to move past the boundary in the first
        // place). To eliminate that false-positive class, re-shuffle
        // against the unclamped `exon_only` bound and only treat the
        // clamp as operative when the unclamped result would have moved
        // strictly past the axis boundary.
        if !matches!(start_axis, boundary::AxisRegion::None) {
            let new_tx_start_0 = new_tx_start.saturating_sub(1);
            let cheap_left = boundaries.left > exon_only.left && new_tx_start_0 == boundaries.left;
            let cheap_right = boundaries.right < exon_only.right && new_tx_end == boundaries.right;
            let direction_could_clamp = match self.config.shuffle_direction {
                ShuffleDirection::FivePrime => cheap_left,
                ShuffleDirection::ThreePrime => cheap_right,
            };
            let (left_clamp_fired, right_clamp_fired) = if direction_could_clamp {
                let (exon_only_start, exon_only_end, _exon_only_edit, _exon_only_warnings) =
                    self.normalize_na_edit(seq, edit, tx_start, tx_end, &exon_only, is_coding)?;
                let exon_only_start_0 = exon_only_start.saturating_sub(1);
                (
                    cheap_left && exon_only_start_0 < boundaries.left,
                    cheap_right && exon_only_end > boundaries.right,
                )
            } else {
                (false, false)
            };
            let (direction_str, clamp_kind) = match self.config.shuffle_direction {
                ShuffleDirection::FivePrime if left_clamp_fired => Some((
                    "5prime",
                    match start_axis {
                        boundary::AxisRegion::Cds => "cds_start",
                        boundary::AxisRegion::ThreeUtr => "3utr",
                        _ => "5utr",
                    },
                )),
                ShuffleDirection::ThreePrime if right_clamp_fired => Some((
                    "3prime",
                    match start_axis {
                        boundary::AxisRegion::Cds => "cds_end",
                        boundary::AxisRegion::FiveUtr => "5utr",
                        _ => "3utr",
                    },
                )),
                _ => None,
            }
            .unwrap_or(("", ""));
            if !direction_str.is_empty() {
                let acc = variant.accession.transcript_accession();
                warnings.push(NormalizationWarning::AxisClampApplied {
                    message: format!(
                        "{}:c.{}: {}-rule shuffle clamped at {} axis boundary; \
                         canonical position was constrained by the CDS\u{2194}UTR sub-axis",
                        acc, variant.loc_edit.location, direction_str, clamp_kind,
                    ),
                    accession: acc,
                    direction: direction_str.to_string(),
                    clamp_kind: clamp_kind.to_string(),
                });
            }
        }

        // Convert back to CDS coordinates
        let new_start = self.tx_to_cds_pos(new_tx_start, cds_start, transcript.cds_end)?;
        let new_end = self.tx_to_cds_pos(new_tx_end, cds_start, transcript.cds_end)?;

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        // Issue #160 + #165 post-canonicalization split. The codon-frame
        // exception applies only to CDS-proper positions, which the
        // helper filters internally via `simple_cds_pos`.
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Cds(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Normalize a transcript variant
    fn normalize_tx(
        &self,
        variant: &TxVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Tx(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins`.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = TxVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return Ok((HV::Tx(new_variant), vec![]));
        }

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Tx(variant.clone()), vec![]));
        }

        // Check for intronic variants or unknown positions
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };

        // Can't normalize variants with unknown (?) offsets - return unchanged
        if has_unknown_offset_tx(start_pos) || has_unknown_offset_tx(end_pos) {
            return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![]));
        }

        // Try to get transcript first - we need it for intronic normalization too.
        // See `normalize_cds` for the rationale behind the dual lookup.
        let accession = variant.accession.transcript_accession();
        let transcript_for_intronic =
            || -> Result<crate::reference::transcript::Transcript, FerroError> {
                self.provider
                    .get_transcript_for_variant(&HV::Tx(variant.clone()))
            };
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            // Can't do full normalization without transcript, but apply minimal notation
            Err(_) => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };

        if start_pos.is_intronic() || end_pos.is_intronic() {
            // Switch to the variant-aware lookup so an NG/NC-parented input
            // gets the build-correct chromosome.
            let transcript = transcript_for_intronic().unwrap_or(transcript);
            // Route intronic tx variants to the intronic normalization path
            if start_pos.is_intronic() && end_pos.is_intronic() {
                return self.normalize_intronic_tx(variant, &transcript, start_pos, end_pos, edit);
            }
            // Variant spans exon-intron boundary - not yet supported for n. coords
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Only normalize positive positions (within transcript)
        // Negative positions are outside the transcript sequence
        if start_pos.base < 1 || end_pos.base < 1 {
            return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![]));
        }

        let tx_start = start_pos.base as u64;
        let tx_end = end_pos.base as u64;

        // Get boundaries
        let boundaries = Boundaries::new(0, transcript.sequence_length());

        // Perform normalization (n. non-coding tx context).
        // Coordinate-only transcripts (no cached bases) fall back to the
        // canonicalize-only path.
        let seq = match transcript.sequence.as_deref() {
            Some(s) => s.as_bytes(),
            None => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };
        let (new_start, new_end, new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, false)?;

        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(TxPos::new(new_start as i64), TxPos::new(new_end as i64)),
                new_edit,
            ),
        };

        // Issue #160 + #165 post-canonicalization split.
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Tx(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Normalize a protein variant
    ///
    /// Protein normalization differs from nucleic acid normalization - there's no
    /// 3'/5' shifting. Instead, we perform formatting standardization:
    ///
    /// 1. **Reference validation**: Check that position amino acids match the
    ///    reference protein sequence (if protein data is available).
    ///
    /// 2. **Redundant sequence removal**: Remove explicit sequences in deletions
    ///    when they match the amino acids at the deletion position.
    ///    Example: `p.Val600delVal` → `p.Val600del`
    ///
    /// 3. **1-letter to 3-letter conversion**: (handled by parser/display)
    fn normalize_protein(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::hgvs::edit::ProteinEdit;
        use crate::hgvs::variant::{LocEdit, ProteinVariant};

        // Validate reference amino acids if provider has protein data
        if self.provider.has_protein_data() {
            self.validate_protein_reference(variant)?;
        }

        // Get the current edit
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Protein(variant.clone()), vec![])),
        };

        // Apply normalization based on edit type
        let normalized_edit = match edit {
            ProteinEdit::Deletion { sequence, count } => {
                // Check for redundant sequence that matches the position
                if let Some(seq) = sequence {
                    if self.is_redundant_protein_deletion_sequence(&variant.loc_edit.location, seq)
                    {
                        // Remove redundant sequence
                        ProteinEdit::Deletion {
                            sequence: None,
                            count: *count,
                        }
                    } else {
                        edit.clone()
                    }
                } else {
                    edit.clone()
                }
            }
            // Other edits pass through unchanged
            _ => edit.clone(),
        };

        // Only create a new variant if the edit changed
        if &normalized_edit != edit {
            let new_variant = ProteinVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| normalized_edit),
                ),
            };
            Ok((HV::Protein(new_variant), vec![]))
        } else {
            Ok((HV::Protein(variant.clone()), vec![]))
        }
    }

    /// Validate that the amino acids in a protein variant match the reference
    ///
    /// Returns an error if there's a mismatch between the variant's stated
    /// amino acid(s) and the actual reference protein sequence.
    fn validate_protein_reference(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
    ) -> Result<(), FerroError> {
        use crate::hgvs::edit::ProteinEdit;

        let accession = variant.accession.transcript_accession();

        // Get start and end positions
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok(()), // Can't validate uncertain positions
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok(()),
        };

        // Get the edit to know what amino acids to validate
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok(()),
        };

        // Only validate edits that specify reference amino acids
        match edit {
            ProteinEdit::Substitution { .. } => {
                // For substitutions, the reference AA comes from the position (start_pos.aa),
                // not from the edit (which may be Xaa placeholder)
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
            }
            ProteinEdit::Deletion { .. } => {
                // Validate start position (from the interval's start AA)
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                // Validate end position if different
                if end_pos.number != start_pos.number {
                    self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
                }
            }
            ProteinEdit::Duplication => {
                // Validate start and end positions
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                if end_pos.number != start_pos.number {
                    self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
                }
            }
            ProteinEdit::Insertion { .. } => {
                // Validate flanking positions
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
            }
            ProteinEdit::Delins { .. } => {
                // Validate start and end
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                if end_pos.number != start_pos.number {
                    self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
                }
            }
            ProteinEdit::Frameshift { .. } => {
                // Validate the frameshift position
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
            }
            ProteinEdit::Extension { .. } => {
                // Extension typically at Ter position - validate
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
            }
            _ => {
                // Identity, Unknown, etc. - no validation needed
            }
        }

        Ok(())
    }

    /// Validate a single protein position against reference
    fn validate_protein_position(
        &self,
        accession: &str,
        position: u64,
        expected_aa: &crate::hgvs::location::AminoAcid,
    ) -> Result<(), FerroError> {
        // Position is 1-based in HGVS, convert to 0-based for sequence access
        // get_protein_sequence uses half-open interval [start, end)
        let start = hgvs_pos_to_index(position) as u64;
        let end = position; // exclusive end

        // Try to get the reference amino acid
        match self.provider.get_protein_sequence(accession, start, end) {
            Ok(ref_seq) => {
                if ref_seq.len() != 1 {
                    return Ok(()); // Unexpected, skip validation
                }

                let ref_aa_char = ref_seq.chars().next().unwrap();
                let expected_char = expected_aa.to_one_letter();

                if ref_aa_char != expected_char {
                    return Err(FerroError::AminoAcidMismatch {
                        accession: accession.to_string(),
                        position,
                        expected: expected_aa.to_string(),
                        found: ref_aa_char.to_string(),
                    });
                }
            }
            Err(_) => {
                // Protein sequence not available, skip validation
            }
        }

        Ok(())
    }

    /// Check if the deletion sequence is redundant (matches the position amino acids)
    ///
    /// A deletion sequence is redundant if it exactly matches the amino acids
    /// specified in the interval. For example:
    /// - `p.Val600delVal` - sequence [Val] matches position 600's Val → redundant
    /// - `p.Lys23_Leu24delLysLeu` - sequence [Lys, Leu] matches positions 23-24 → redundant
    fn is_redundant_protein_deletion_sequence(
        &self,
        interval: &crate::hgvs::interval::ProtInterval,
        sequence: &crate::hgvs::edit::AminoAcidSeq,
    ) -> bool {
        // Get the start and end positions
        let start_pos = match interval.start.inner() {
            Some(pos) => pos,
            None => return false,
        };
        let end_pos = match interval.end.inner() {
            Some(pos) => pos,
            None => return false,
        };

        // Calculate expected sequence length from interval
        let interval_len = if end_pos.number >= start_pos.number {
            (end_pos.number - start_pos.number + 1) as usize
        } else {
            return false;
        };

        // Check if sequence length matches
        if sequence.len() != interval_len {
            return false;
        }

        // For a point deletion (single AA), check if the sequence matches the position AA
        if interval_len == 1 {
            return sequence.0.len() == 1 && sequence.0[0] == start_pos.aa;
        }

        // For a range deletion, check first and last AAs
        // The sequence should be [start_aa, ..., end_aa]
        if let (Some(first), Some(last)) = (sequence.0.first(), sequence.0.last()) {
            return *first == start_pos.aa && *last == end_pos.aa;
        }

        false
    }

    /// Normalize an RNA variant
    ///
    /// RNA variants (r.) are similar to transcript variants (n.) and undergo
    /// the same 3'/5' shifting normalization for indels. The main difference
    /// is that RNA uses lowercase nucleotides in HGVS notation.
    fn normalize_rna(
        &self,
        variant: &crate::hgvs::variant::RnaVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::variant::{LocEdit, RnaVariant};

        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins`.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = RnaVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return Ok((HV::Rna(new_variant), vec![]));
        }

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Rna(variant.clone()), vec![]));
        }

        // Check for intronic variants or unknown positions
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        if start_pos.is_intronic() || end_pos.is_intronic() {
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Try to get transcript (RNA uses the same accession as mRNA transcripts)
        let accession = variant.accession.transcript_accession();
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            Err(_) => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // Convert RNA positions to transcript-1 positions, deciding per
        // endpoint. UTR (`r.*N`/`r.-N`) and non-positive bases need a CDS
        // to translate; non-UTR positive bases map 1:1 to transcript-1
        // indices and work without a CDS (mock providers used in tests
        // often omit cds_start/end). Choosing per endpoint keeps mixed
        // intervals like `r.50_*1del` from rerouting the positive end
        // through `rna_to_tx_pos` — issue #163 follow-up.
        //
        // Axis convention (pinned by `tests/issue_291_rna_axis_convention.rs`,
        // closes #291): `r.` positive non-UTR bases are
        // **transcript-1-relative**, NOT CDS-relative. `r.10` against a
        // transcript with `cds_start = 100` maps to tx index 10, not tx 109.
        // This is consistent across `fetch_ref_for_canonical_split`,
        // `simple_rna_pos`, and `normalize_na_edit`. Only `r.*N`/`r.-N`
        // (and base 0) translate through `cds_start`/`cds_end` via
        // `rna_to_tx_pos`.
        let cds_info = transcript.cds_start.zip(transcript.cds_end);
        let map_in = |pos: &crate::hgvs::location::RnaPos| -> Option<u64> {
            if pos.utr3 || pos.base < 1 {
                let (cds_start, cds_end) = cds_info?;
                self.rna_to_tx_pos(pos, cds_start, Some(cds_end)).ok()
            } else {
                Some(pos.base as u64)
            }
        };
        let tx_start = match map_in(start_pos) {
            Some(v) => v,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };
        let tx_end = match map_in(end_pos) {
            Some(v) => v,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // Get boundaries (entire transcript span; r. has no exon-level
        // junction restriction beyond the transcript ends).
        let boundaries = Boundaries::new(0, transcript.sequence_length());

        // Perform normalization (RNA context: codon-frame gate does not apply;
        // r. is not in the spec's accepted reference types for repeats).
        // Coordinate-only transcripts fall back to the canonicalize-only path.
        let seq = match transcript.sequence.as_deref() {
            Some(s) => s.as_bytes(),
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };
        let (new_tx_start, new_tx_end, new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, false)?;

        // Convert each normalized tx position back independently, restoring
        // UTR notation when the position falls outside the CDS. This catches
        // both the original issue #163 case (UTR input shuffling within the
        // UTR) and a positive-base input that shifts past `cds_end` during
        // normalization. Without `cds_info` we keep the simple base-1 mapping.
        //
        // For positions that stay inside the CDS-proper window
        // (`cds_start <= pos <= cds_end`), the tx index is emitted directly
        // as `r.{pos}` — i.e. the transcript-1-relative axis is preserved on
        // the way out, matching the convention enforced by `map_in` above.
        // See `tests/issue_291_rna_axis_convention.rs` for the pin.
        use crate::hgvs::location::RnaPos;
        let map_out = |pos: u64| -> Result<RnaPos, FerroError> {
            if let Some((cds_start, cds_end)) = cds_info {
                if pos < cds_start || pos > cds_end {
                    return self.tx_to_rna_pos(pos, cds_start, Some(cds_end));
                }
            }
            Ok(RnaPos::new(pos as i64))
        };
        let new_start = map_out(new_tx_start)?;
        let new_end = map_out(new_tx_end)?;

        let new_variant = RnaVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(RnaInterval::new(new_start, new_end), new_edit),
        };

        // Issue #160 + #165 post-canonicalization split (T/U-equivalent
        // comparison for the rev-comp scan and per-position emissions).
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Rna(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Convert an RNA position to a transcript-1 position.
    ///
    /// Mirrors `cds_to_tx_pos`. `r.*N` maps to `cds_end + N`, `r.-N` to
    /// `cds_start + N` (HGVS skips the `0` gap).
    fn rna_to_tx_pos(
        &self,
        pos: &crate::hgvs::location::RnaPos,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<u64, FerroError> {
        if pos.utr3 {
            let end = cds_end.ok_or_else(|| FerroError::ConversionError {
                msg: "No CDS end".to_string(),
            })?;
            let base = u64::try_from(pos.base).map_err(|_| FerroError::ConversionError {
                msg: format!("Negative base {} in 3' UTR position", pos.base),
            })?;
            Ok(end + base)
        } else if pos.base < 0 {
            let tx_pos = cds_start as i64 + pos.base;
            u64::try_from(tx_pos).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "RNA position r.{} maps before transcript start (cds_start={})",
                    pos.base, cds_start
                ),
            })
        } else if pos.base == 0 {
            Ok(cds_start.saturating_sub(1))
        } else {
            Ok(cds_start + pos.base as u64 - 1)
        }
    }

    /// Convert a transcript-1 position back to an RNA position, restoring
    /// the appropriate region (`r.*N` for 3'UTR, `r.-N` for 5'UTR).
    fn tx_to_rna_pos(
        &self,
        pos: u64,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<crate::hgvs::location::RnaPos, FerroError> {
        use crate::hgvs::location::RnaPos;
        let end = cds_end.ok_or_else(|| FerroError::ConversionError {
            msg: "No CDS end".to_string(),
        })?;
        if pos < cds_start {
            Ok(RnaPos {
                base: pos as i64 - cds_start as i64,
                offset: None,
                utr3: false,
            })
        } else if pos > end {
            Ok(RnaPos {
                base: (pos - end) as i64,
                offset: None,
                utr3: true,
            })
        } else {
            Ok(RnaPos {
                base: (pos - cds_start + 1) as i64,
                offset: None,
                utr3: false,
            })
        }
    }

    /// Normalize a mitochondrial variant
    ///
    /// Mirrors `normalize_genome` for non-origin-crossing variants: fetch
    /// a sequence window around the variant, run `normalize_na_edit` with
    /// `is_coding=false` (mito is genomic-style and not subject to the
    /// codon-frame restriction — `repeated.md` line 21 restricts the
    /// codon-frame gate exclusively to `c.` descriptions), then map
    /// positions back.
    ///
    /// Origin-crossing (wraparound) `del`/`delins` variants are rejected
    /// at parse time by `parse_genome_interval`'s inverted-range check;
    /// wraparound `dup`/`ins`/`inv` are exempt from that check and reach
    /// this function. Without provider data they fall through to
    /// `canonicalize_mt_variant`; with provider data the window fetch
    /// errors (start > end is invalid for a linear slice), again
    /// dropping to the fallback. F1 / #129 will introduce circular-aware
    /// semantics in a follow-up — see `tests/mito_circular_audit.rs` for
    /// the pinned behavior preserved by this PR.
    fn normalize_mt(
        &self,
        variant: &crate::hgvs::variant::MtVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions.
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Mt(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins` up front. Pure-syntax;
        // no reference data needed. Returns immediately (no inv-split
        // pass), matching `normalize_genome`: callers that want
        // inv-split on the rewritten delins can re-normalize.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = MtVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return Ok((HV::Mt(new_variant), vec![]));
        }

        // Only normalize indels; substitutions / identity / repeat-with-
        // count pass through unchanged. Mirrors `normalize_genome`.
        if !needs_normalization(edit) {
            return Ok((HV::Mt(variant.clone()), vec![]));
        }

        // Fallback for variants we cannot remap through the window-based
        // pipeline (unknown position, decorated position, or no provider
        // data). Runs minimal-notation cleanup, then still applies
        // `apply_canonical_split` so issue #160 inv-split and issue #165
        // sub-only decomposition remain in force — those run on a narrow
        // fetch (`fetch_ref_for_canonical_split`) that can succeed even
        // when the wider shuffle window does not.
        let mt_fallback = |v: &crate::hgvs::variant::MtVariant| {
            let canonical = self.canonicalize_mt_variant(v);
            let (split, warnings) = self.apply_canonical_split(HV::Mt(canonical));
            (wrap_allele_if_split(split), warnings)
        };

        let accession = variant.accession.transcript_accession();
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok(mt_fallback(variant)),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok(mt_fallback(variant)),
        };

        // Decorated genome positions (offset / pter / qter / cen) cannot
        // be losslessly remapped through base-only window normalization —
        // remapping via `pos.base` and rebuilding with `GenomePos::new`
        // would silently drop the decoration. Fall back to minimal-
        // notation cleanup for these.
        if start_pos.offset.is_some()
            || end_pos.offset.is_some()
            || start_pos.is_special()
            || end_pos.is_special()
        {
            return Ok(mt_fallback(variant));
        }

        let start = start_pos.base;
        let end = end_pos.base;

        // Window-based fetch around the variant. Non-origin-crossing
        // variants take this path exactly like genomic; wraparound
        // `dup`/`ins`/`inv` (where `start > end`) drop through to the
        // canonicalize-only fallback below via `get_sequence` failure.
        let window_start = start.saturating_sub(self.config.window_size);
        let seq_result = self.provider.get_sequence(
            &accession,
            window_start,
            end.saturating_add(self.config.window_size),
        );

        let ref_seq = match seq_result {
            Ok(s) => s,
            // No reference data → fall back to minimal-notation cleanup.
            Err(_) => return Ok(mt_fallback(variant)),
        };

        let rel_start = start - window_start;
        let rel_end = end - window_start;

        // Mitochondrial reference is plus-strand and not subject to the
        // codon-frame `unit_len % 3 == 0` restriction (the mito genome
        // has no canonical "CDS" exemption boundary in the same sense as
        // nuclear `c.`; the spec's mito chapter doesn't carry the
        // codon-frame clause), so pass `is_coding=false`.
        let (new_rel_start, new_rel_end, new_edit, mut warnings) = self.normalize_na_edit(
            ref_seq.as_bytes(),
            edit,
            rel_start,
            rel_end,
            &Boundaries::new(0, ref_seq.len() as u64),
            false,
        )?;

        let new_start = new_rel_start + window_start;
        let new_end = new_rel_end + window_start;

        let new_variant = MtVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(new_start), GenomePos::new(new_end)),
                new_edit,
            ),
        };

        // Issue #160 inv-split post-pass (mirrors normalize_genome).
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Mt(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Apply minimal notation to an mt variant without full normalization.
    /// Mirrors `canonicalize_genome_variant` — used as a fallback when
    /// reference data is unavailable.
    fn canonicalize_mt_variant(
        &self,
        variant: &crate::hgvs::variant::MtVariant,
    ) -> crate::hgvs::variant::MtVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        MtVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Normalize an intronic CDS variant
    ///
    /// This converts the intronic position to genomic coordinates, normalizes
    /// in genomic space, and converts back to CDS intronic notation.
    fn normalize_intronic_cds(
        &self,
        variant: &CdsVariant,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &CdsPos,
        end_pos: &CdsPos,
        edit: &NaEdit,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::convert::CoordinateMapper;

        // Check if we have genomic data available
        if !self.provider.has_genomic_data() {
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Get the chromosome for this transcript. Issue #332: include the
        // parent accession and full variant Display in the error so the
        // remaining failure mode (transcript not present on any genome build)
        // is diagnosable without re-running with extra logging.
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!(
                        "Transcript {} has no chromosome mapping for intronic \
                         normalization (parent={}, variant={}); the cdot data has \
                         no genomic alignment for this transcript on any known \
                         genome build",
                        transcript.id,
                        variant
                            .accession
                            .genomic_context
                            .as_deref()
                            .map(|a| a.full())
                            .unwrap_or_else(|| "<none>".to_string()),
                        variant,
                    ),
                })?;

        // Create coordinate mapper
        let mapper = CoordinateMapper::new(transcript);

        // Convert CDS intronic positions to genomic
        let g_start = mapper.cds_to_genomic_with_intron(start_pos)?;
        let g_end = mapper.cds_to_genomic_with_intron(end_pos)?;

        // On minus strand, genomic coords may be reversed relative to coding order.
        // Track whether we swap so we can restore coding order after normalization.
        let swapped = g_start > g_end;
        let (g_start, g_end) = if swapped {
            (g_end, g_start)
        } else {
            (g_start, g_end)
        };

        // Get a window of genomic sequence around the variant for normalization
        // Use the same window size as for exonic normalization
        let window = self.config.window_size;
        let seq_start = g_start.saturating_sub(window);
        let seq_end = g_end.saturating_add(window);

        // Fetch genomic sequence
        let genomic_seq = self
            .provider
            .get_genomic_sequence(chromosome, seq_start, seq_end)?;

        // Calculate the variant position relative to the fetched sequence
        let rel_start = (g_start - seq_start) + 1; // 1-based
        let rel_end = (g_end - seq_start) + 1;

        // Define boundaries within the intron
        // For intronic variants, we can shift within the intron but not into exons
        // Find the intron boundaries
        // Use exon-aware CDS-to-tx mapping to account for cdot's gap positions
        let tx_pos = mapper.cds_to_tx(start_pos)?;
        let tx_start = u64::try_from(tx_pos.base).map_err(|_| FerroError::ConversionError {
            msg: format!(
                "Negative transcript position {} during intronic normalization",
                tx_pos.base
            ),
        })?;

        let intron = transcript
            .find_intron_at_tx_boundary(tx_start, start_pos.offset.unwrap_or(0))
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find intron for normalization".to_string(),
            })?;

        // Get intron boundaries in genomic coordinates
        let (intron_g_start, intron_g_end) = match (intron.genomic_start, intron.genomic_end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                return Err(FerroError::ConversionError {
                    msg: "Intron has no genomic coordinates".to_string(),
                })
            }
        };

        // Calculate relative intron boundaries
        let intron_rel_start = intron_g_start.saturating_sub(seq_start) + 1;
        let intron_rel_end = intron_g_end.saturating_sub(seq_start) + 1;
        let boundaries = Boundaries::new(intron_rel_start, intron_rel_end);

        // On minus-strand transcripts the genomic-strand sequence is the
        // reverse complement of the transcript view, but the variant's
        // edit alt is in transcript view. Running `normalize_na_edit` on
        // the genomic-strand bytes therefore defeats every rule that
        // compares the alt against the local reference window. Flip the
        // sequence and the relative positions / boundaries to transcript
        // view here, run normalization, then map the result positions
        // back to the genomic frame. (Issue #98.)
        let (work_seq, work_rel_start, work_rel_end, work_boundaries) = flip_intronic_for_strand(
            transcript.strand,
            &genomic_seq,
            rel_start,
            rel_end,
            &boundaries,
        );

        // Perform normalization in transcript-view space (CDS intronic context).
        // HGVS spec (repeated.md): the codon-frame restriction
        // (`unit_len % 3 == 0` for repeat notation in `c.` context)
        // applies only to bases inside the CDS proper. Introns are
        // exempt:
        //   > This restriction only applies to the coding sequence,
        //   > which does not include the introns or the UTR sequence.
        // Pass `is_coding=false` so an intronic homopolymer dup/del
        // can emit `[N±k]` repeat notation instead of falling back to
        // the gated `ins<literal>` / plain `del` forms.
        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            false,
        )?;

        // Map the result positions back to the genomic-strand frame
        let (new_rel_start, new_rel_end) = unflip_intronic_positions(
            transcript.strand,
            work_seq.len() as u64,
            work_new_rel_start,
            work_new_rel_end,
        );

        // Convert the normalized genomic position back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert back to CDS intronic notation
        let new_start = mapper.genomic_to_cds_intronic(new_g_start)?;
        let new_end = mapper.genomic_to_cds_intronic(new_g_end)?;

        // Restore coding order if positions were swapped for genomic processing
        let (new_start, new_end) = if swapped {
            (new_end, new_start)
        } else {
            (new_start, new_end)
        };

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        Ok((HV::Cds(new_variant), warnings))
    }

    /// Normalize an intronic transcript (n.) variant
    ///
    /// This mirrors `normalize_intronic_cds()` but works with TxPos instead of CdsPos.
    /// Converts to genomic coordinates, normalizes in genomic space, and converts back.
    fn normalize_intronic_tx(
        &self,
        variant: &TxVariant,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &TxPos,
        end_pos: &TxPos,
        edit: &NaEdit,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::convert::CoordinateMapper;

        // Check if we have genomic data available
        if !self.provider.has_genomic_data() {
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Get the chromosome for this transcript. Issue #332: include the
        // parent accession and full variant Display in the error so the
        // remaining failure mode (transcript not present on any genome build)
        // is diagnosable without re-running with extra logging.
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!(
                        "Transcript {} has no chromosome mapping for intronic \
                         normalization (parent={}, variant={}); the cdot data has \
                         no genomic alignment for this transcript on any known \
                         genome build",
                        transcript.id,
                        variant
                            .accession
                            .genomic_context
                            .as_deref()
                            .map(|a| a.full())
                            .unwrap_or_else(|| "<none>".to_string()),
                        variant,
                    ),
                })?;

        // Create coordinate mapper
        let mapper = CoordinateMapper::new(transcript);

        // Convert tx intronic positions to genomic
        let g_start = mapper.tx_to_genomic_with_intron(start_pos)?;
        let g_end = mapper.tx_to_genomic_with_intron(end_pos)?;

        // On minus strand, genomic coords may be reversed relative to coding order.
        // Track whether we swap so we can restore coding order after normalization.
        let swapped = g_start > g_end;
        let (g_start, g_end) = if swapped {
            (g_end, g_start)
        } else {
            (g_start, g_end)
        };

        // Get a window of genomic sequence around the variant
        let window = self.config.window_size;
        let seq_start = g_start.saturating_sub(window);
        let seq_end = g_end.saturating_add(window);

        // Fetch genomic sequence
        let genomic_seq = self
            .provider
            .get_genomic_sequence(chromosome, seq_start, seq_end)?;

        // Calculate the variant position relative to the fetched sequence
        let rel_start = (g_start - seq_start) + 1; // 1-based
        let rel_end = (g_end - seq_start) + 1;

        // Find the intron boundaries for normalization limits
        let tx_start = u64::try_from(start_pos.base).map_err(|_| FerroError::ConversionError {
            msg: format!(
                "Negative transcript position {} during intronic normalization",
                start_pos.base
            ),
        })?;

        let intron = transcript
            .find_intron_at_tx_boundary(tx_start, start_pos.offset.unwrap_or(0))
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find intron for normalization".to_string(),
            })?;

        // Get intron boundaries in genomic coordinates
        let (intron_g_start, intron_g_end) = match (intron.genomic_start, intron.genomic_end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                return Err(FerroError::ConversionError {
                    msg: "Intron has no genomic coordinates".to_string(),
                })
            }
        };

        // Calculate relative intron boundaries
        let intron_rel_start = intron_g_start.saturating_sub(seq_start) + 1;
        let intron_rel_end = intron_g_end.saturating_sub(seq_start) + 1;
        let boundaries = Boundaries::new(intron_rel_start, intron_rel_end);

        // See `normalize_intronic_cds`: same orientation fix for #98.
        let (work_seq, work_rel_start, work_rel_end, work_boundaries) = flip_intronic_for_strand(
            transcript.strand,
            &genomic_seq,
            rel_start,
            rel_end,
            &boundaries,
        );

        // Perform normalization in transcript-view space (n. non-coding intronic context).
        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            false,
        )?;

        // Map the result positions back to the genomic-strand frame
        let (new_rel_start, new_rel_end) = unflip_intronic_positions(
            transcript.strand,
            work_seq.len() as u64,
            work_new_rel_start,
            work_new_rel_end,
        );

        // Convert the normalized genomic position back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert back to transcript intronic notation
        let new_start = mapper.genomic_to_tx_with_intron(new_g_start)?;
        let new_end = mapper.genomic_to_tx_with_intron(new_g_end)?;

        // Restore coding order if positions were swapped for genomic processing
        let (new_start, new_end) = if swapped {
            (new_end, new_start)
        } else {
            (new_start, new_end)
        };

        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        Ok((HV::Tx(new_variant), warnings))
    }

    /// Normalize a CDS variant that spans an exon-intron boundary
    ///
    /// This handles cases like:
    /// - c.914_918+3del (exonic start, intronic end)
    /// - c.194-64_233del (intronic start, exonic end)
    ///
    /// Strategy: Convert to genomic coordinates, normalize there, convert back.
    fn normalize_boundary_spanning_cds(
        &self,
        variant: &CdsVariant,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &CdsPos,
        end_pos: &CdsPos,
        edit: &NaEdit,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::convert::CoordinateMapper;

        // Require genomic data for boundary spanning normalization
        if !self.provider.has_genomic_data() {
            return Err(FerroError::ExonIntronBoundary {
                exon: 0,
                variant: format!("{}", variant),
            });
        }

        // Issue #332: same improved error shape as the intronic paths.
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!(
                        "Transcript {} has no chromosome for boundary \
                         normalization (parent={}, variant={}); the cdot data \
                         has no genomic alignment for this transcript on any \
                         known genome build",
                        transcript.id,
                        variant
                            .accession
                            .genomic_context
                            .as_deref()
                            .map(|a| a.full())
                            .unwrap_or_else(|| "<none>".to_string()),
                        variant,
                    ),
                })?;

        let mapper = CoordinateMapper::new(transcript);

        // Convert both positions to genomic
        // For exonic positions, use standard conversion
        // For intronic positions, use intronic conversion
        let g_start = self.cds_pos_to_genomic(&mapper, start_pos)?;
        let g_end = self.cds_pos_to_genomic(&mapper, end_pos)?;

        // On minus strand, genomic coords may be reversed relative to coding order.
        // Track whether we swap so we can restore coding order after normalization.
        let swapped = g_start > g_end;
        let (g_start, g_end) = if swapped {
            (g_end, g_start)
        } else {
            (g_start, g_end)
        };

        // Fetch genomic sequence with window for normalization
        let window = self.config.window_size;
        let seq_start = g_start.saturating_sub(window);
        let seq_end = g_end.saturating_add(window);
        let genomic_seq = self
            .provider
            .get_genomic_sequence(chromosome, seq_start, seq_end)?;

        // Calculate relative positions (1-based)
        let rel_start = (g_start - seq_start) + 1;
        let rel_end = (g_end - seq_start) + 1;

        // Determine normalization boundaries
        // For boundary-spanning variants, use the union of exon and intron boundaries
        let boundaries =
            self.get_boundary_spanning_limits(transcript, &mapper, start_pos, end_pos, seq_start)?;

        // On minus-strand transcripts the genomic-strand window is the
        // reverse complement of the transcript view, but the variant's edit
        // alt is in transcript view. Running `normalize_na_edit` on raw
        // genomic bytes therefore canonicalizes against the wrong alphabet
        // (and the codon-frame repeat gate inspects ref context here too).
        // Mirror the intronic flow: flip into transcript view before
        // normalization, then unflip the result positions back to the
        // genomic frame. (CDS boundary-spanning context.)
        let (work_seq, work_rel_start, work_rel_end, work_boundaries) = flip_intronic_for_strand(
            transcript.strand,
            &genomic_seq,
            rel_start,
            rel_end,
            &boundaries,
        );

        // HGVS spec (repeated.md): the codon-frame restriction applies
        // only to bases inside the CDS proper. Boundary-spanning
        // variants cross an exon/intron boundary, so their footprint
        // is not entirely within coding sequence — pass
        // `is_coding=false` to match the intronic exemption. (A
        // hypothetical purely-exonic-span variant won't enter this
        // function; the exonic CDS path in `normalize_cds` makes its
        // own UTR/CDS-aware choice.)
        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            false,
        )?;

        let (new_rel_start, new_rel_end) = unflip_intronic_positions(
            transcript.strand,
            work_seq.len() as u64,
            work_new_rel_start,
            work_new_rel_end,
        );

        // Convert back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert genomic back to CDS
        // The result might be:
        // - Still boundary-spanning
        // - Fully exonic (if shifted into exon)
        // - Fully intronic (if shifted into intron)
        let new_start = mapper.genomic_to_cds_intronic(new_g_start)?;
        let new_end = mapper.genomic_to_cds_intronic(new_g_end)?;

        // Restore coding order if positions were swapped for genomic processing
        let (new_start, new_end) = if swapped {
            (new_end, new_start)
        } else {
            (new_start, new_end)
        };

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        Ok((HV::Cds(new_variant), warnings))
    }

    /// Convert a CDS position (exonic or intronic) to genomic coordinate
    fn cds_pos_to_genomic(
        &self,
        mapper: &crate::convert::CoordinateMapper,
        pos: &CdsPos,
    ) -> Result<u64, FerroError> {
        if pos.is_intronic() {
            mapper.cds_to_genomic_with_intron(pos)
        } else {
            mapper
                .cds_to_genomic(pos)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("CDS position {} not in exons", pos.base),
                })
        }
    }

    /// Calculate normalization boundaries for boundary-spanning variants
    ///
    /// Returns the genomic region within which we can shift the variant,
    /// which is the union of the exon and intron containing the variant ends.
    fn get_boundary_spanning_limits(
        &self,
        transcript: &crate::reference::transcript::Transcript,
        mapper: &crate::convert::CoordinateMapper,
        start_pos: &CdsPos,
        end_pos: &CdsPos,
        seq_start: u64,
    ) -> Result<Boundaries, FerroError> {
        // Find the genomic extent of the region we can shift within
        // This is the union of:
        // - The exon containing the exonic position
        // - The intron containing the intronic position

        // Identify which position is exonic and which is intronic
        let (exonic_pos, intronic_pos) = if start_pos.is_intronic() {
            (end_pos, start_pos)
        } else {
            (start_pos, end_pos)
        };

        // Get exon boundaries in genomic coords
        let tx_pos = mapper.cds_to_tx(exonic_pos)?;
        let tx_pos_base = u64::try_from(tx_pos.base).map_err(|_| FerroError::ConversionError {
            msg: format!(
                "Negative transcript position {} during boundary normalization",
                tx_pos.base
            ),
        })?;
        let exon = transcript
            .exon_at(tx_pos_base)
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find exon for boundary normalization".to_string(),
            })?;

        // Get intron boundaries in genomic coords
        let tx_boundary = mapper.cds_to_tx(intronic_pos)?;
        let tx_boundary_base =
            u64::try_from(tx_boundary.base).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "Negative transcript position {} during boundary normalization",
                    tx_boundary.base
                ),
            })?;
        let offset = intronic_pos.offset.unwrap_or(0);
        let intron = transcript
            .find_intron_at_tx_boundary(tx_boundary_base, offset)
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find intron for boundary normalization".to_string(),
            })?;

        // Get genomic coordinates from exon and intron
        let (exon_g_start, exon_g_end): (u64, u64) = match (exon.genomic_start, exon.genomic_end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                return Err(FerroError::ConversionError {
                    msg: "Exon has no genomic coordinates".to_string(),
                })
            }
        };

        let (intron_g_start, intron_g_end): (u64, u64) =
            match (intron.genomic_start, intron.genomic_end) {
                (Some(s), Some(e)) => (s, e),
                _ => {
                    return Err(FerroError::ConversionError {
                        msg: "Intron has no genomic coordinates".to_string(),
                    })
                }
            };

        // Union of exon and intron genomic coordinates
        let combined_start = exon_g_start.min(intron_g_start);
        let combined_end = exon_g_end.max(intron_g_end);

        // Convert to relative positions (1-based within our fetched sequence)
        let rel_start = combined_start.saturating_sub(seq_start) + 1;
        let rel_end = combined_end.saturating_sub(seq_start) + 1;

        Ok(Boundaries::new(rel_start, rel_end))
    }

    /// Core normalization for nucleic acid edits
    ///
    /// Returns (new_start, new_end, new_edit, warnings)
    fn normalize_na_edit(
        &self,
        ref_seq: &[u8],
        edit: &NaEdit,
        start: u64,
        end: u64,
        boundaries: &Boundaries,
        is_coding: bool,
    ) -> Result<(u64, u64, NaEdit, Vec<NormalizationWarning>), FerroError> {
        let mut warnings = Vec::new();

        // Validate reference allele before normalization
        let validation = validate::validate_reference(edit, ref_seq, start, end);
        if !validation.valid {
            // `corrected` is honest about whether the canonical Display
            // drops the user-stated bases. See issue #280.
            //
            // True for `Deletion` / `Duplication` (canonicalize_edit
            // strips `sequence`) and conditionally for `Inversion`
            // (the Inversion arm in `normalize_na_edit` emits
            // `sequence: None` when `shorten_inversion` rewrites the
            // span, and otherwise still drops the stated bases in the
            // canonical path). `Substitution` does not reach here at
            // all — `needs_normalization` returns `false` for real
            // substitutions and only the degenerate `ref == alt` case
            // routes through this function, where the validator's
            // single-base check is satisfied by construction.
            //
            // False for `Repeat` / `MultiRepeat` consistency mismatches
            // (issues #214 / #279): the per-unit declaration is part of
            // the user's form and the normalizer passes the description
            // through verbatim.
            //
            // TODO/Note: revisit if `delins` ever gets a stated-deleted
            // validator — today `validate_reference`'s `NaEdit::Delins`
            // arm returns `ok()` unconditionally (see also the Delins
            // arm in `canonicalize_edit`, which strips `deleted` /
            // `deleted_length`), so the branch is unreachable for
            // Delins.
            let corrected = !matches!(edit, NaEdit::Repeat { .. } | NaEdit::MultiRepeat { .. });
            warnings.push(NormalizationWarning::RefSeqMismatch {
                message: validation.warning.unwrap_or_default(),
                stated_ref: validation.stated_ref.unwrap_or_default(),
                actual_ref: validation.actual_ref.unwrap_or_default(),
                position: format!("{}-{}", start, end),
                corrected,
            });
        }

        // Substitution with ref == alt is identity (e.g. c.100A>A → c.100=).
        // This is the SNV companion to the same-base delins → identity rule;
        // the rewrite is purely syntactic on the edit's stated bases, so it
        // applies across coordinate systems and runs before shuffling.
        if let NaEdit::Substitution {
            reference,
            alternative,
        } = edit
        {
            if reference == alternative {
                return Ok((start, end, NaEdit::position_identity(), warnings));
            }
        }

        // Get the alternate sequence for the edit
        let alt_seq = match edit {
            NaEdit::Deletion { .. } => vec![],
            NaEdit::Insertion { sequence } => {
                // Only shuffle if we have a literal sequence
                if let Some(bases) = sequence.bases() {
                    bases.iter().map(|b| b.to_u8()).collect()
                } else {
                    // Cannot shuffle non-literal insertions (counts, ranges, etc.)
                    return Ok((start, end, edit.clone(), warnings.clone()));
                }
            }
            NaEdit::Duplication { .. } => {
                // Always read duplicated bytes from the reference,
                // regardless of any user-stated bases on the input
                // (`dup<base>` shapes). The canonical HGVS Display
                // drops the stated bases anyway (see the dup output
                // arm in `get_canonical_form`); using the stated
                // bases here when they don't match the reference
                // (the parser accepts them in lenient mode with a
                // `RefSeqMismatch` warning) caused the shuffle to
                // mis-shift and broke single-pass idempotency:
                // `c.10dupA` against a `CCCCC` homopolymer
                // canonicalized to `c.10dup` on pass 1 (stated-ref
                // stripped without shifting) and only shifted to
                // `c.14dup` on pass 2 (because the now-clean form
                // reads bytes from reference and the shuffle fires
                // correctly). Reading from reference up-front
                // collapses both pass 1 and pass 2 behavior into a
                // single pass. Issue #219.
                let s = hgvs_pos_to_index(start);
                let e = end as usize;
                if e <= ref_seq.len() {
                    ref_seq[s..e].to_vec()
                } else {
                    vec![]
                }
            }
            NaEdit::Delins { sequence, .. } => {
                use crate::hgvs::edit::InsertedSequence;

                // HGVS spec (issue #81 A3): a delins with an empty inserted
                // sequence is semantically a deletion and must be rendered as
                // `del`. Rewrite up-front so the result picks up del 3'-shift
                // and validation in the Deletion arm.
                if matches!(sequence, InsertedSequence::Empty) {
                    let del = NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    };
                    return self
                        .normalize_na_edit(ref_seq, &del, start, end, boundaries, is_coding);
                }

                // HGVS spec: delins should NOT be 3' shifted like del/dup/ins,
                // but the edit-type priority (sub > del > inv > dup > ins) means
                // we may need to rewrite it as a higher-priority form: identity
                // (insert == ref), substitution (1->1 ref!=alt), or duplication.
                if let InsertedSequence::Literal(seq) = sequence {
                    use crate::hgvs::edit::{Base, Sequence};
                    use rules::DelinsCanonical;
                    let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();
                    let start_idx = hgvs_pos_to_index(start);
                    let end_idx = end as usize;

                    // Reconstruct an InsertedSequence from a Vec<u8> produced by
                    // shared-affix trimming. The bytes round-trip through `Base`
                    // because they originated from a typed `Sequence` (the input
                    // `seq` above), so `from_char` cannot fail; expect-on-None
                    // makes the invariant explicit if a future refactor breaks
                    // the pipeline.
                    let bytes_to_inserted_seq = |bytes: &[u8]| -> InsertedSequence {
                        let bases: Vec<Base> = bytes
                            .iter()
                            .map(|b| {
                                Base::from_char(*b as char).expect(
                                    "trimmed delins byte must be a valid IUPAC base \
                                     because the input sequence was already a typed Sequence",
                                )
                            })
                            .collect();
                        InsertedSequence::Literal(Sequence::new(bases))
                    };

                    match rules::canonicalize_delins(ref_seq, start_idx, end_idx, &seq_bytes) {
                        DelinsCanonical::Identity => {
                            // c.10delinsG where ref[10]=G  ->  c.10=
                            return Ok((start, end, NaEdit::position_identity(), warnings.clone()));
                        }
                        DelinsCanonical::Substitution {
                            position,
                            reference,
                            alternative,
                        } => {
                            // g.1000delinsA where ref[1000]=G  ->  g.1000G>A.
                            // After shared-affix trimming `position` is a
                            // 0-indexed offset into ref_seq, not necessarily
                            // the input `start`.
                            let pos = index_to_hgvs_pos(position);
                            return Ok((
                                pos,
                                pos,
                                NaEdit::Substitution {
                                    reference,
                                    alternative,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::Deletion { start: s0, end: e0 } => {
                            // c.2_5delinsAT (ref ACGT) -> c.3_4del. Range fields
                            // are the trimmed half-open 0-indexed interval.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Deletion {
                                    sequence: None,
                                    length: None,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::Insertion {
                            after_index,
                            sequence: ins_bytes,
                        } => {
                            // c.2_4delinsACGT (ref ACT) -> c.3_4insG. `after_index`
                            // is the 0-indexed position of the base AFTER the
                            // insertion, which is also the 1-based HGVS position
                            // of the base BEFORE — so HGVS X = after_index,
                            // Y = after_index + 1.
                            //
                            // Recurse into `normalize_na_edit` with the new
                            // Insertion so the full ins pipeline (3'/5' shuffle
                            // + `insertion_to_duplication` + `insertion_to_repeat`)
                            // runs. Without recursion, a delins-derived
                            // insertion that duplicates a nearby reference tract
                            // (e.g. biocommons `g.X delinsCTTTCTT` where
                            // ref[X+1..X+6]=TTTCTT) would skip the
                            // ins→dup recognizer and emit the long `insTTTCTT`
                            // form instead of the canonical `dup`. Closes-after:
                            // #356.
                            let new_edit = NaEdit::Insertion {
                                sequence: bytes_to_inserted_seq(&ins_bytes),
                            };
                            // Preserve warnings collected for the original
                            // delins (e.g. RefSeqMismatch in strict mode):
                            // merge them into the recursive call's warnings
                            // rather than dropping them by tail-returning.
                            let (new_start, new_end, new_edit, mut child_warnings) = self
                                .normalize_na_edit(
                                    ref_seq,
                                    &new_edit,
                                    after_index as u64,
                                    (after_index + 1) as u64,
                                    boundaries,
                                    is_coding,
                                )?;
                            let mut merged = warnings.clone();
                            merged.append(&mut child_warnings);
                            return Ok((new_start, new_end, new_edit, merged));
                        }
                        DelinsCanonical::Inversion { start: s0, end: e0 } => {
                            // A2 (#81): g.100_102delinsTAG where ref=CTA  ->  g.100_102inv.
                            // Position interval is already shortened. e0 is the
                            // exclusive 0-based end; the HGVS 1-based inclusive
                            // end takes the same numeric value.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Inversion {
                                    sequence: None,
                                    length: None,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::Duplication { start: s0, end: e0 } => {
                            // c.5delinsGG where ref[5]=G  ->  c.5dup. Duplication
                            // is detected before trimming, so the range matches
                            // the input.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Duplication {
                                    sequence: None,
                                    length: None,
                                    uncertain_extent: None,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::KeepAsDelins {
                            start: s0,
                            end: e0,
                            sequence: trimmed_bytes,
                        } => {
                            // Either no trimming was possible (range == input)
                            // or trimming reduced the delins to a smaller delins
                            // that still doesn't fit a higher-priority form.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Delins {
                                    sequence: bytes_to_inserted_seq(&trimmed_bytes),
                                    deleted: None,
                                    deleted_length: None,
                                },
                                warnings.clone(),
                            ));
                        }
                    }
                }
                // Non-literal insert (Count, Range, Reference, PositionRange,
                // Complex, …): cannot trim or classify without the actual
                // bases, but we still strip an explicit deleted sequence /
                // length per the same spec rule that the Literal arm above
                // applies (`delins.md`: the recommendation is to omit the
                // explicit deleted bases). Closes #338's WITH-provider gap
                // surfaced in code review.
                return Ok((start, end, canonicalize_edit(edit), warnings.clone()));
            }
            NaEdit::Inversion { .. } => {
                // Apply the complementary-outer-bases shortening rule. After
                // shortening, the inversion's interval may no longer match the
                // input's explicit `sequence`/`length` (if any), so emit
                // minimal notation.
                let start_idx = hgvs_pos_to_index(start); // Convert 1-based to 0-based
                let end_idx = end as usize; // end is exclusive (0-based)

                if let Some((new_s, new_e)) = rules::shorten_inversion(ref_seq, start_idx, end_idx)
                {
                    return Ok((
                        index_to_hgvs_pos(new_s),
                        new_e as u64,
                        NaEdit::Inversion {
                            sequence: None,
                            length: None,
                        },
                        warnings,
                    ));
                } else {
                    // Inversion reduced to identity. Use the canonical
                    // position-only Identity (matches the Delins identity arm
                    // above), so both inversion-collapse paths emit the same
                    // shape.
                    return Ok((start, end, NaEdit::position_identity(), warnings));
                }
            }
            NaEdit::Repeat {
                sequence,
                count,
                additional_counts,
                trailing,
            } => {
                use crate::hgvs::edit::RepeatCount;

                // Only normalize exact counts with a sequence
                let Some(seq) = sequence else {
                    return Ok((start, end, edit.clone(), warnings.clone()));
                };

                // Range, UncertainRange, MinUncertain, MaxUncertain, Unknown:
                // no concrete count to 3'-shift against, so pass through unchanged.
                let RepeatCount::Exact(specified_count) = count else {
                    return Ok((start, end, edit.clone(), warnings.clone()));
                };

                // Skip if there are additional counts (genotype notation)
                if !additional_counts.is_empty() || trailing.is_some() {
                    return Ok((start, end, edit.clone(), warnings.clone()));
                }

                // Get the repeat unit as bytes
                let repeat_unit: Vec<u8> = seq.bases().iter().map(|b| b.to_u8()).collect();
                let pos_idx = hgvs_pos_to_index(start); // Convert 1-based to 0-based

                // Normalize the repeat
                match rules::normalize_repeat(
                    ref_seq,
                    pos_idx,
                    &repeat_unit,
                    *specified_count,
                    is_coding,
                ) {
                    rules::RepeatNormResult::Deletion {
                        start: del_start,
                        end: del_end,
                    } => {
                        // Minimal notation - no explicit length
                        let del_edit = NaEdit::Deletion {
                            sequence: None,
                            length: None,
                        };
                        return Ok((del_start, del_end, del_edit, warnings));
                    }
                    rules::RepeatNormResult::Duplication {
                        start: dup_start,
                        end: dup_end,
                        sequence: _dup_seq,
                    } => {
                        // Minimal notation - no explicit sequence or length
                        let dup_edit = NaEdit::Duplication {
                            sequence: None,
                            length: None,
                            uncertain_extent: None,
                        };
                        return Ok((dup_start, dup_end, dup_edit, warnings));
                    }
                    rules::RepeatNormResult::Insertion {
                        start: ins_start,
                        end: ins_end,
                        sequence: ins_seq,
                    } => {
                        // Codon-frame gate routed an expansion to ins literal form
                        // (e.g., c.1741_1742insTATATATA per spec).
                        use crate::hgvs::edit::{Base, InsertedSequence, Sequence};
                        let bases: Vec<Base> = ins_seq
                            .iter()
                            .filter_map(|&b| Base::from_char(b as char))
                            .collect();
                        if bases.len() == ins_seq.len() {
                            let ins_edit = NaEdit::Insertion {
                                sequence: InsertedSequence::Literal(Sequence::new(bases)),
                            };
                            return Ok((ins_start, ins_end, ins_edit, warnings));
                        }
                        // Defensive fallback: rule layer returned a base byte
                        // that doesn't fit the Base alphabet (e.g. N). Don't
                        // emit a truncated insertion — keep the original edit
                        // and positions so downstream invariants hold.
                        return Ok((start, end, edit.clone(), warnings));
                    }
                    rules::RepeatNormResult::Repeat {
                        start: rep_start,
                        end: rep_end,
                        sequence: rep_seq,
                        count: rep_count,
                    } => {
                        use crate::hgvs::edit::{Base, RepeatCount, Sequence};
                        let bases: Vec<Base> = rep_seq
                            .iter()
                            .filter_map(|&b| Base::from_char(b as char))
                            .collect();
                        if bases.len() == rep_seq.len() {
                            let rep_edit = NaEdit::Repeat {
                                sequence: Some(Sequence::new(bases)),
                                count: RepeatCount::Exact(rep_count),
                                additional_counts: vec![],
                                trailing: None,
                            };
                            return Ok((rep_start, rep_end, rep_edit, warnings));
                        }
                        // Defensive fallback: rule layer returned a repeat
                        // unit byte that doesn't fit the Base alphabet (e.g.
                        // a gap or non-IUPAC byte from the reference). Don't
                        // emit a truncated repeat sequence — keep the
                        // original edit and positions.
                        return Ok((start, end, edit.clone(), warnings));
                    }
                    rules::RepeatNormResult::Unchanged => {
                        return Ok((start, end, edit.clone(), warnings.clone()));
                    }
                }
            }
            _ => return Ok((start, end, edit.clone(), warnings.clone())), // Other edits don't need shuffling
        };

        // Perform shuffle
        // For insertions, the HGVS interval X_Y (where Y = X+1) represents flanking positions.
        // We need to adjust the end coordinate so shuffle checks the correct reference position.
        // For c.445_446insA: start=634, end=635 (1-based tx coords)
        // We want shuffle to check ref_seq[634-1] = ref_seq[633] for first flanking
        // and ref_seq[635-1] = ref_seq[634] for second flanking (the position to check for 3' shift)
        //
        // For insertions, we need to adjust the start passed to shuffle so that the alt_idx
        // calculation starts at 0 (not 1). The shuffle's alt_idx formula is:
        //   alt_idx = (new_end - start) % alt_seq.len()
        // For insertions, new_end starts at shuffle_end (which is end - 1 for insertions).
        // If we pass start_idx directly, alt_idx = (end-1) - start_idx = 1 (wrong, should be 0).
        // By passing start_idx + 1, we get alt_idx = (end-1) - (start_idx+1) = 0 (correct).
        let shuffle_end = match edit {
            NaEdit::Insertion { .. } => end.saturating_sub(1), // Use second flanking position
            _ => end,
        };
        let start_idx = hgvs_pos_to_index(start); // Convert 1-based to 0-based
        let shuffle_start = match edit {
            NaEdit::Insertion { .. } => start_idx as u64 + 1, // Adjust so alt_idx starts at 0
            _ => start_idx as u64,
        };
        let result = shuffle(
            ref_seq,
            &alt_seq,
            shuffle_start,
            shuffle_end, // Adjusted for insertions
            boundaries,
            self.config.shuffle_direction,
        );

        // Convert back to 1-based HGVS positions
        // For insertions, we adjusted the start for shuffle, now adjust back
        let shuffle_result_start = match edit {
            NaEdit::Insertion { .. } => result.start.saturating_sub(1), // Adjust back
            _ => result.start,
        };
        let new_start = index_to_hgvs_pos(shuffle_result_start as usize);
        // For insertions, we adjusted the end for shuffle, now restore the HGVS X_(X+1) format
        let new_end = match edit {
            NaEdit::Insertion { .. } => index_to_hgvs_pos(result.end as usize), // Restore second flanking position
            _ => result.end,
        };

        // Determine the canonical form for the edit
        // HGVS rules:
        // - Deletions ALWAYS stay as deletions (just shift 3')
        // - Insertions become duplications if single-base matches adjacent
        // - Multi-base insertions/dups in homopolymer become repeat notation
        let (final_start, final_end, new_edit) = match edit {
            NaEdit::Insertion { sequence } => {
                use crate::hgvs::edit::{InsertedSequence, RepeatCount, Sequence};

                if let InsertedSequence::Literal(seq) = sequence {
                    let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();

                    // Check for repeat notation first (multi-base homopolymer insertion)
                    // Use the ORIGINAL position (start), not shuffled position (result.start)
                    // because repeat notation refers to the reference tract position
                    if seq_bytes.len() > 1 {
                        let original_pos_idx = hgvs_pos_to_index(start) as u64; // 0-based original position
                        if let Some((_first, count, rep_start, rep_end, unit_bytes)) =
                            rules::insertion_to_repeat(
                                ref_seq,
                                original_pos_idx,
                                &seq_bytes,
                                is_coding,
                            )
                        {
                            use crate::hgvs::edit::Base;
                            let bases: Vec<Base> = unit_bytes
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            if bases.len() == unit_bytes.len() {
                                let repeat_seq = Sequence::new(bases);
                                let repeat_edit = NaEdit::Repeat {
                                    sequence: Some(repeat_seq),
                                    count: RepeatCount::Exact(count),
                                    additional_counts: vec![],
                                    trailing: None,
                                };
                                return Ok((rep_start, rep_end, repeat_edit, warnings));
                            }
                        }
                    }

                    // Resolve insertion → duplication canonicalization. Three candidate
                    // dup positions compete; we pick by the rules below in order.
                    //
                    // (a) Tract-aligned dup via `insertion_to_duplication` (uses the
                    //     ORIGINAL insertion point and finds the maximal tandem run
                    //     under any cyclic rotation of the alt). When the tract has
                    //     `ref_count >= 2` we prefer this regardless of how far
                    //     shuffle walked: the multi-copy tract has a meaningful phase
                    //     that the spec-canonical form preserves (issue #132).
                    //
                    // (b) Post-shuffle simple dup via `insertion_is_duplication`. When
                    //     shuffle walked past a single-copy tract via partial-match
                    //     extension (e.g. TGATC abutting TGAAG — first three bases
                    //     match but the fourth does not), the post-shuffle position is
                    //     more 3' than (a)'s tract-aligned position and is the canonical
                    //     answer per the 3' rule (issue #180).
                    //
                    // (c) Single-copy tract fallback (`insertion_to_duplication` with
                    //     `ref_count == 1`). Hit when shuffle stalled before completing
                    //     one alt rotation (so (b) doesn't find a dup at the post-
                    //     shuffle position) but the alt does match an adjacent ref unit
                    //     at the ORIGINAL position. Example: ins AACA abutting AACA.
                    //
                    // If none match, fall through to ins (possibly rotated).
                    let original_pos_idx = hgvs_pos_to_index(start) as u64;
                    let ins_to_dup = rules::insertion_to_duplication(
                        ref_seq,
                        original_pos_idx,
                        &seq_bytes,
                        self.config.shuffle_direction,
                    );

                    // Codon-frame gate (repeated.md): in c., if the alt is
                    // >=2 copies of a non-codon-aligned unit, the spec
                    // mandates ins<literal>, not dup. The smallest-unit
                    // length is rotation-invariant, so we can compute it
                    // once from `seq_bytes` and apply the same gate at
                    // every dup-emission site below. In practice the gate
                    // never fires when `ins_to_dup` is `Some` (that helper
                    // only returns for single-unit alts, where
                    // `smallest_unit.len() == seq_bytes.len()`), but we
                    // guard the (a) fast path and (c) single-copy fallback
                    // anyway so the spec rule is explicit at each site and
                    // survives future changes to `insertion_to_duplication`.
                    let smallest_unit = rules::smallest_repeat_unit(&seq_bytes);
                    let codon_blocks_dup = is_coding
                        && smallest_unit.len() < seq_bytes.len()
                        && !smallest_unit.len().is_multiple_of(3);

                    if !codon_blocks_dup {
                        if let Some(rules::InsToDupResult {
                            start: dup_start,
                            end: dup_end,
                            ref_count,
                            ..
                        }) = ins_to_dup.as_ref()
                        {
                            if *ref_count >= 2 {
                                return Ok((
                                    *dup_start,
                                    *dup_end,
                                    NaEdit::Duplication {
                                        sequence: None,
                                        length: None,
                                        uncertain_extent: None,
                                    },
                                    warnings,
                                ));
                            }
                        }
                    }

                    // Check for simple duplication (single-base or matching adjacent)
                    // When shifting an insertion through a repeat region, the effective sequence
                    // rotates. For example, shifting "GGC" by 1 position gives "GCG".
                    let shift_amount =
                        (result.start as usize).saturating_sub(shuffle_start as usize);
                    let rotation = shift_amount % seq_bytes.len();
                    let rotated_seq: Vec<u8> = if rotation > 0 {
                        seq_bytes[rotation..]
                            .iter()
                            .chain(seq_bytes[..rotation].iter())
                            .copied()
                            .collect()
                    } else {
                        seq_bytes.clone()
                    };

                    if !codon_blocks_dup
                        && rules::insertion_is_duplication(ref_seq, result.start, &rotated_seq)
                    {
                        // For duplication, use minimal notation without explicit sequence
                        // Position: for c.X_(X+1)ins that duplicates preceding sequence,
                        // the result is c.(X-len+1)_Xdup
                        let dup_len = rotated_seq.len() as u64;
                        let (dup_start, dup_end) = if dup_len == 1 {
                            (new_start, new_start) // Single position for single-base dup
                        } else {
                            // Multi-base dup: the duplicated region is BEFORE the insertion point
                            (new_start - dup_len + 1, new_start)
                        };
                        (
                            dup_start,
                            dup_end,
                            NaEdit::Duplication {
                                sequence: None, // Minimal notation - no explicit sequence
                                length: None,
                                uncertain_extent: None,
                            },
                        )
                    } else if let Some(rules::InsToDupResult {
                        start: dup_start,
                        end: dup_end,
                        ..
                    }) = ins_to_dup.as_ref().filter(|_| !codon_blocks_dup)
                    {
                        // (c) Single-copy tract fallback. Reached when (a) declined
                        // because `ref_count < 2` and (b) declined because the post-
                        // shuffle rotated alt doesn't match adjacent reference. The
                        // alt is a (possibly rotated) tandem unit abutting a single-
                        // copy ref tract at the original insertion point — emit the
                        // dup over that tract.
                        (
                            *dup_start,
                            *dup_end,
                            NaEdit::Duplication {
                                sequence: None,
                                length: None,
                                uncertain_extent: None,
                            },
                        )
                    } else {
                        // Output the rotated sequence for shifted insertions
                        if rotation > 0 {
                            use crate::hgvs::edit::{Base, InsertedSequence, Sequence};
                            let rotated_bases: Vec<Base> = rotated_seq
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            // Mirror the gated-ins guard used by the
                            // RepeatNormResult::Insertion / GatedInsertion
                            // branches: if any byte fell outside the Base
                            // alphabet, refuse to emit a truncated `ins`
                            // and fall back to the original edit so
                            // downstream invariants hold.
                            if rotated_bases.len() == rotated_seq.len() {
                                let new_sequence =
                                    InsertedSequence::Literal(Sequence::new(rotated_bases));
                                (
                                    new_start,
                                    new_end,
                                    NaEdit::Insertion {
                                        sequence: new_sequence,
                                    },
                                )
                            } else {
                                (new_start, new_end, edit.clone())
                            }
                        } else {
                            (new_start, new_end, edit.clone())
                        }
                    }
                } else {
                    (new_start, new_end, edit.clone())
                }
            }
            NaEdit::Duplication { .. } => {
                use crate::hgvs::edit::{Base, RepeatCount, Sequence};

                // Check if duplication should become repeat notation
                // Use the shuffled positions (result.start, result.end) which are 0-based
                // This applies to both single-base dups in homopolymers and multi-base tandem dups
                if let Some(dup_result) =
                    rules::duplication_to_repeat(ref_seq, result.start, result.end, is_coding)
                {
                    match dup_result {
                        rules::DupToRepeatResult::Homopolymer {
                            base,
                            count,
                            start: rep_start,
                            end: rep_end,
                        } => {
                            if let Some(base_enum) = Base::from_char(base as char) {
                                let repeat_seq = Sequence::new(vec![base_enum]);
                                let repeat_edit = NaEdit::Repeat {
                                    sequence: Some(repeat_seq),
                                    count: RepeatCount::Exact(count),
                                    additional_counts: vec![],
                                    trailing: None,
                                };
                                return Ok((rep_start, rep_end, repeat_edit, warnings));
                            }
                        }
                        rules::DupToRepeatResult::TandemRepeat {
                            unit,
                            count,
                            start: rep_start,
                            end: rep_end,
                        } => {
                            let bases: Vec<Base> = unit
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            if bases.len() == unit.len() {
                                let repeat_seq = Sequence::new(bases);
                                let repeat_edit = NaEdit::Repeat {
                                    sequence: Some(repeat_seq),
                                    count: RepeatCount::Exact(count),
                                    additional_counts: vec![],
                                    trailing: None,
                                };
                                return Ok((rep_start, rep_end, repeat_edit, warnings));
                            }
                        }
                        rules::DupToRepeatResult::GatedInsertion {
                            start: ins_start,
                            end: ins_end,
                            sequence: ins_seq,
                        } => {
                            // Codon-frame gate routed a multi-copy dup to ins
                            // literal form per HGVS spec.
                            use crate::hgvs::edit::InsertedSequence;
                            let bases: Vec<Base> = ins_seq
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            if bases.len() == ins_seq.len() {
                                let ins_edit = NaEdit::Insertion {
                                    sequence: InsertedSequence::Literal(Sequence::new(bases)),
                                };
                                return Ok((ins_start, ins_end, ins_edit, warnings));
                            }
                            // Defensive fallback: rule layer returned a base
                            // byte that doesn't fit the Base alphabet (e.g.
                            // N). Fall through to the generic dup minimal-
                            // notation path below rather than emitting a
                            // truncated insertion.
                        }
                    }
                }
                // Keep as duplication but strip explicit sequence (minimal notation)
                (
                    new_start,
                    new_end,
                    NaEdit::Duplication {
                        sequence: None,
                        length: None,
                        uncertain_extent: None,
                    },
                )
            }
            // Deletions: post-shift, check for B2 canonical-form rule
            // (deletion of >=2 tandem-repeat units → unit[N-k]); otherwise
            // strip explicit length for minimal `del` notation. The
            // collect-into-Option short-circuits if any byte in the unit isn't
            // a valid `Base` (e.g. `N`), in which case we fall through to del.
            //
            // B2 is defined for a *post-3'-shift* deletion (the shuffle phase-
            // alignment lemma justifies emitting `unit[N-k]` without rotation).
            // Under FivePrime shuffle, applying it would re-anchor the
            // 5'-normalized deletion to the canonical tract position, defeating
            // the user's choice of direction — so gate it on ThreePrime.
            NaEdit::Deletion { .. } => {
                use crate::hgvs::edit::{Base, RepeatCount, Sequence};
                if self.config.shuffle_direction == ShuffleDirection::ThreePrime {
                    if let Some(rep) = rules::deletion_to_repeat(
                        ref_seq,
                        result.start as usize,
                        result.end as usize,
                        is_coding,
                    ) {
                        let bases: Option<Vec<Base>> = rep
                            .unit
                            .iter()
                            .map(|&b| Base::from_char(b as char))
                            .collect();
                        if let Some(bases) = bases {
                            let repeat_edit = NaEdit::Repeat {
                                sequence: Some(Sequence::new(bases)),
                                count: RepeatCount::Exact(rep.count),
                                additional_counts: vec![],
                                trailing: None,
                            };
                            return Ok((rep.start, rep.end, repeat_edit, warnings));
                        }
                    }
                }
                (
                    new_start,
                    new_end,
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                )
            }
            // All other edit types stay unchanged
            _ => (new_start, new_end, edit.clone()),
        };

        Ok((final_start, final_end, new_edit, warnings))
    }

    /// Convert CDS position to transcript position
    fn cds_to_tx_pos(
        &self,
        pos: &CdsPos,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<u64, FerroError> {
        if pos.utr3 {
            let end = cds_end.ok_or_else(|| FerroError::ConversionError {
                msg: "No CDS end".to_string(),
            })?;
            let base = u64::try_from(pos.base).map_err(|_| FerroError::ConversionError {
                msg: format!("Negative base {} in 3' UTR position", pos.base),
            })?;
            Ok(end + base)
        } else if pos.base < 0 {
            // 5'UTR: HGVS numbering skips c.0 (c.-1 is the base immediately
            // upstream of c.1), so c.-N maps to tx position cds_start - N.
            // Issue #97 — the previous formula `cds_start + base - 1`
            // double-counted the gap and emitted the wrong tx position.
            let tx_pos = cds_start as i64 + pos.base;
            u64::try_from(tx_pos).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "CDS position c.{} maps before transcript start (cds_start={})",
                    pos.base, cds_start
                ),
            })
        } else if pos.base == 0 {
            // c.0 is not a valid HGVS position, but historical inputs
            // can land here. Preserve the legacy mapping (treat as the
            // last 5'UTR base, equivalent to c.-1) rather than failing.
            Ok(cds_start.saturating_sub(1))
        } else {
            Ok(cds_start + pos.base as u64 - 1)
        }
    }

    /// Convert transcript position to CDS position
    fn tx_to_cds_pos(
        &self,
        pos: u64,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<CdsPos, FerroError> {
        let end = cds_end.ok_or_else(|| FerroError::ConversionError {
            msg: "No CDS end".to_string(),
        })?;

        if pos < cds_start {
            // 5'UTR: HGVS numbering skips c.0, so a tx position one
            // base 5' of cds_start is c.-1 (not c.0). Inverse of the
            // forward formula `tx = cds_start + base` for negative
            // base: `base = tx - cds_start`. Issue #97 — the previous
            // formula `tx - cds_start + 1` would emit base = 0 for
            // tx = cds_start - 1, rendered by `CdsPos::Display` as
            // `c.?` (`CDS_BASE_UNKNOWN`).
            Ok(CdsPos {
                base: pos as i64 - cds_start as i64,
                offset: None,
                utr3: false,
            })
        } else if pos > end {
            Ok(CdsPos {
                base: (pos - end) as i64,
                offset: None,
                utr3: true,
            })
        } else {
            Ok(CdsPos {
                base: (pos - cds_start + 1) as i64,
                offset: None,
                utr3: false,
            })
        }
    }

    /// Apply minimal notation to a CDS variant without full normalization.
    ///
    /// This is used when we can't do full normalization (e.g., missing transcript)
    /// but still want to apply minimal HGVS notation rules.
    fn canonicalize_cds_variant(&self, variant: &CdsVariant) -> CdsVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        // Only canonicalize if the edit has redundant information
        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Apply minimal notation to a genome variant without full normalization.
    fn canonicalize_genome_variant(&self, variant: &GenomeVariant) -> GenomeVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        GenomeVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Apply minimal notation to a transcript variant without full normalization.
    fn canonicalize_tx_variant(&self, variant: &TxVariant) -> TxVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Post-canonicalization split for a single normalized variant.
    /// Coord-system-agnostic: handles `g.`, `m.`, `c.` (CDS-proper positions
    /// only), `n.`, and `r.`. Fetches the per-coord-system reference window
    /// internally, calls `decompose_delins`, and rebuilds N variants when
    /// the decomposition fires. Returns `vec![variant]` if the variant
    /// doesn't decompose (non-Delins, complex location, no provider data,
    /// nothing to split out, **or the provider's ref window doesn't match
    /// the HGVS coordinate span** — e.g. cdot alignment gaps shorten the
    /// returned byte slice, see #339).
    ///
    /// Implements two spec-priority rules from `general.md:56`
    /// (substitution > deletion > inversion > duplication > insertion):
    /// - Inversion priority: a delins whose span contains a rev-comp
    ///   sub-span splits into `[…; inv; …]` (issue #160).
    /// - Substitution priority: a delins whose post-trim span contains
    ///   two or more independent single-base mismatches separated by at
    ///   least one unchanged nucleotide splits into separate substitutions
    ///   (issue #165 / item A10). The narrow codon-frame exception
    ///   (`general.md:35-38`) is preserved by `build_split_variants`,
    ///   which re-groups `[Sub; Identity; Sub]` triplets whose endpoints
    ///   share a codon when the variant is in CDS.
    ///
    /// Position math: `decompose_delins` returns 0-indexed offsets into
    /// the fetched `ref_bytes` slice. `ref_bytes[0]` corresponds to the
    /// variant's HGVS start position, so absolute HGVS pos = `hgvs_start +
    /// offset`.
    ///
    /// RNA `r.` variants have `U` bases in the alt while transcript ref bytes
    /// are `T`. Both slices are normalized to `T` before comparison so the
    /// rev-comp scan works uniformly; the emitted `Substitution` sub-edits
    /// preserve the original alt `Base` (which may be `Base::U`).
    fn apply_canonical_split(
        &self,
        variant: HgvsVariant,
    ) -> (Vec<HgvsVariant>, Vec<NormalizationWarning>) {
        let Some((hgvs_start, hgvs_end, alt_bytes, ref_bytes)) =
            self.fetch_ref_for_canonical_split(&variant)
        else {
            return (vec![variant], vec![]);
        };
        // When the provider returns fewer (or more) bytes than the HGVS
        // interval span, the canonical-split decomposition would walk past
        // the end of `ref_bytes` (or under-walk and miss interior
        // identities). This happens in practice when a cdot exon
        // alignment collapses gaps that the HGVS span counts as
        // positions (e.g. biocommons NG_032871.1:g.32476_53457delins…
        // returns 15,539 bytes for a 20,982-bp span). Pre-fix this
        // fired a `debug_assert_eq!` that panicked debug builds and
        // would have produced out-of-bounds reads in release. Now bail
        // out gracefully with a `CanonicalSplitSkipped` warning so
        // callers can flag the input for human review. Closes #339,
        // closes-after #354.
        let n = ref_bytes.len();
        let expected_span = (hgvs_end - hgvs_start + 1) as usize;
        if n != expected_span {
            let accession = variant_accession_string(&variant);
            let warning = NormalizationWarning::CanonicalSplitSkipped {
                message: format!(
                    "{}: variant span {}..{} ({} bp) exceeds provider's reference window \
                     ({} bp). Per HGVS spec refseq.md \u{00A7}43, the variant must be \
                     entirely encompassed by the reference. Strict mode rejects.",
                    accession, hgvs_start, hgvs_end, expected_span, n,
                ),
                accession,
                hgvs_start,
                hgvs_end,
                expected_span,
                actual_bytes: n,
            };
            return (vec![variant], vec![warning]);
        }
        let ref_norm = normalize_t_u(&ref_bytes);
        let alt_norm = normalize_t_u(&alt_bytes);
        let Some(subedits) = rules::decompose_delins(&ref_norm, 0, n, &alt_norm) else {
            return (vec![variant], vec![]);
        };
        // Substitution sub-edits inherit `alt_norm` bytes (T-form) from
        // `decompose_delins`, but the user's literal alt may have been U
        // (r. inputs). Re-derive the substitution `alternative` from the
        // pre-normalized `alt_bytes` so r. variants render `g>u` instead of
        // a silently coerced `g>t`. The position field is a 0-indexed offset
        // into the same window passed to `decompose_delins`, so it indexes
        // alt_bytes directly.
        let subedits = subedits
            .into_iter()
            .map(|se| match se {
                rules::DelinsSubedit::Substitution {
                    position,
                    reference,
                    alternative,
                } => {
                    let alt = crate::hgvs::edit::Base::from_char(alt_bytes[position] as char)
                        .unwrap_or(alternative);
                    rules::DelinsSubedit::Substitution {
                        position,
                        reference,
                        alternative: alt,
                    }
                }
                other => other,
            })
            .collect();
        // The codon-frame exception (`general.md:35-38`) applies to any
        // variant on the CDS-relative axis: `c.` (CDS proper) and `r.`
        // (RNA CDS-relative; issue #275 item 1).
        // `fetch_ref_for_canonical_split` already filters to CDS-proper
        // positions via `simple_cds_pos` / `simple_rna_pos`, so the
        // discriminant check below is sufficient. The exception fires
        // inside `build_split_variants` for every embedded
        // `[Sub; Identity; Sub]` triplet whose endpoints share a codon.
        // Caveat: the `r.` branch of `fetch_ref_for_canonical_split`
        // does not translate the RNA axis through `cds_start`, so for
        // transcripts with `cds_start > 1` the ref window can be
        // mis-aligned. The codon-frame split path inherits that latent
        // bug. See #291 — the fix is to apply `rna_to_tx_pos` (or
        // equivalent) before reading bytes.
        let codon_frame_aware = matches!(variant, HgvsVariant::Cds(_) | HgvsVariant::Rna(_));
        (
            build_split_variants(&variant, subedits, hgvs_start, codon_frame_aware),
            vec![],
        )
    }

    /// Per-coord-system extraction of `(hgvs_start, hgvs_end, alt_bytes,
    /// ref_bytes)` for the post-canonicalization split. Returns `None`
    /// when the variant is not a single-Delins at simple positions, or
    /// when the provider can't supply the ref window.
    ///
    /// The `ref_bytes` slice is sized exactly to the variant's HGVS interval
    /// (`hgvs_end - hgvs_start + 1` bytes), with `ref_bytes[0]` aligned to
    /// HGVS pos `hgvs_start`. This invariant lets the caller use a uniform
    /// `hgvs_pos = hgvs_start + offset` formula regardless of coord system.
    ///
    /// Note: this helper reads the transcript *sequence* (FASTA-derived,
    /// build-invariant) and not the `chromosome` field; the bare
    /// `provider.get_sequence(&transcript_accession, …)` is therefore
    /// sufficient. The build-aware lookup (`get_transcript_for_variant`,
    /// see #332) is reserved for paths that actually consume `chromosome`
    /// — the intronic and boundary-spanning normalization branches.
    fn fetch_ref_for_canonical_split(
        &self,
        variant: &HgvsVariant,
    ) -> Option<(u64, u64, Vec<u8>, Vec<u8>)> {
        let (hgvs_start, hgvs_end, alt) = extract_simple_delins(variant)?;
        let ref_bytes = match variant {
            HgvsVariant::Genome(g) => self
                .provider
                // get_sequence is 0-based half-open: [hgvs_start - 1, hgvs_end).
                .get_sequence(
                    &g.accession.transcript_accession(),
                    hgvs_start - 1,
                    hgvs_end,
                )
                .ok()?
                .into_bytes(),
            HgvsVariant::Mt(m) => self
                .provider
                .get_sequence(
                    &m.accession.transcript_accession(),
                    hgvs_start - 1,
                    hgvs_end,
                )
                .ok()?
                .into_bytes(),
            HgvsVariant::Cds(c) => {
                // CDS pos N → 1-based tx pos = cds_start + N - 1.
                // 0-based tx slice = [cds_start + N - 2, cds_start + end - 1).
                let tx = self
                    .provider
                    .get_transcript(&c.accession.transcript_accession())
                    .ok()?;
                let cds_start = tx.cds_start?;
                let s = cds_start.checked_add(hgvs_start)?.checked_sub(2)? as usize;
                let e = cds_start.checked_add(hgvs_end)?.checked_sub(1)? as usize;
                let bytes = tx.sequence.as_deref()?.as_bytes();
                if e > bytes.len() || s >= e {
                    return None;
                }
                bytes[s..e].to_vec()
            }
            HgvsVariant::Tx(t) => {
                let tx = self
                    .provider
                    .get_transcript(&t.accession.transcript_accession())
                    .ok()?;
                let s = (hgvs_start - 1) as usize;
                let e = hgvs_end as usize;
                let bytes = tx.sequence.as_deref()?.as_bytes();
                if e > bytes.len() || s >= e {
                    return None;
                }
                bytes[s..e].to_vec()
            }
            HgvsVariant::Rna(r) => {
                // `r.` positive non-UTR bases use the transcript-1-relative
                // axis in this codebase (NOT CDS-relative): `r.N` for `N > 0`
                // and non-UTR maps directly to tx index `N`. This matches the
                // convention used by `simple_rna_pos`, `normalize_rna::map_in`/
                // `map_out`, and `normalize_na_edit` for the r. arm.
                // `simple_rna_pos` filters intronic positions, 3'-UTR
                // positions (the `r.*N` form), and non-positive bases. It
                // does NOT filter 5'-UTR positions because under the tx-1
                // axis convention those have positive bases below
                // `cds_start` and are themselves valid tx indices — the
                // `hgvs_start - 1` slice against the full transcript
                // sequence is therefore correct for any position reaching
                // this arm (including tx-1 5'-UTR).
                // Pinned by `tests/issue_291_rna_axis_convention.rs` (closes #291).
                let tx = self
                    .provider
                    .get_transcript(&r.accession.transcript_accession())
                    .ok()?;
                let s = (hgvs_start - 1) as usize;
                let e = hgvs_end as usize;
                let bytes = tx.sequence.as_deref()?.as_bytes();
                if e > bytes.len() || s >= e {
                    return None;
                }
                bytes[s..e].to_vec()
            }
            _ => return None,
        };
        Some((hgvs_start, hgvs_end, alt, ref_bytes))
    }

    /// Issue #160 + #165 post-merge canonicalization for a single
    /// variant. Used by the cis-allele merge path; `normalize_allele`
    /// applies this per merged variant. Conservatively returns
    /// `vec![v]` for variants the helper can't process.
    ///
    /// Three spec rules are folded together by re-running normalization
    /// on the merged variant:
    /// - Full-span canonicalization (identity / dup / sub / del / ins /
    ///   full-span inv with outer-pair shortening) handled by
    ///   `canonicalize_delins` inside `normalize_na_edit`.
    /// - Sub-span inv decomposition (the issue #160 case) handled by
    ///   `apply_canonical_split` wired into each per-coord-system
    ///   `normalize_*`.
    /// - Sub-only decomposition for delins containing interior identities
    ///   (issue #165 / item A10), with the spec's codon-frame exception
    ///   (`general.md:35-38`) preserved inside `build_split_variants`.
    ///
    /// If the result is an `HgvsVariant::Allele` (the split fired and
    /// produced multiple variants), unwrap its inner variants so they
    /// flatten into the outer cis-allele list rather than nesting.
    fn canonical_split_for_variant(&self, v: HgvsVariant) -> Vec<HgvsVariant> {
        if !matches!(
            v,
            HgvsVariant::Genome(_)
                | HgvsVariant::Mt(_)
                | HgvsVariant::Cds(_)
                | HgvsVariant::Tx(_)
                | HgvsVariant::Rna(_)
        ) {
            return vec![v];
        }
        if extract_simple_delins(&v).is_none() {
            return vec![v];
        }
        match self.normalize_with_warnings(&v) {
            Ok(r) => match r.result {
                HgvsVariant::Allele(a) => a.variants,
                other => vec![other],
            },
            Err(_) => vec![v],
        }
    }
}

/// Flip a fetched intronic genomic-strand window into transcript-view
/// orientation when the host transcript is on the minus strand. Returns
/// the input unchanged on plus strand. The relative positions and the
/// shuffle boundaries are flipped so they index into the returned
/// sequence consistently. Companion to [`unflip_intronic_positions`].
fn flip_intronic_for_strand(
    strand: Strand,
    genomic_seq: &str,
    rel_start: u64,
    rel_end: u64,
    boundaries: &Boundaries,
) -> (String, u64, u64, Boundaries) {
    if strand != Strand::Minus {
        return (
            genomic_seq.to_string(),
            rel_start,
            rel_end,
            boundaries.clone(),
        );
    }
    let seq_len = genomic_seq.len() as u64;
    let rc = crate::sequence::reverse_complement(genomic_seq);
    let new_rel_start = seq_len - rel_end + 1;
    let new_rel_end = seq_len - rel_start + 1;
    let new_boundaries = Boundaries::new(
        seq_len - boundaries.right + 1,
        seq_len - boundaries.left + 1,
    );
    (rc, new_rel_start, new_rel_end, new_boundaries)
}

/// Inverse of [`flip_intronic_for_strand`] for the result positions
/// emitted by `normalize_na_edit`. On plus strand returns the input
/// unchanged; on minus strand maps from transcript-view back to the
/// genomic-strand frame.
fn unflip_intronic_positions(
    strand: Strand,
    seq_len: u64,
    rel_start: u64,
    rel_end: u64,
) -> (u64, u64) {
    if strand == Strand::Minus {
        (seq_len - rel_end + 1, seq_len - rel_start + 1)
    } else {
        (rel_start, rel_end)
    }
}

// =============================================================================
// Issue #160 + #165: delins post-canonicalization split helpers
// =============================================================================
//
// After `normalize_na_edit` (or `merge_consecutive_edits` for cis alleles)
// produces a Delins variant, the resulting span may be expressible in a
// higher-priority form under `general.md:56` (sub > del > inv > dup > ins).
// Two cases fire here:
// - Inversion sub-span: the delins span contains a rev-comp sub-region —
//   split into `[…; inv; …]` (issue #160).
// - Independent substitutions: the delins span contains two or more
//   single-base mismatches separated by at least one unchanged nucleotide
//   — split each into its own sub variant (issue #165 / tracking issue
//   #81 item A10). The codon-frame exception (`general.md:35-38`) is
//   preserved when applicable (see `build_split_variants`).
//
// The split is implemented as a post-pass over an already-built variant. It
// fetches a reference window via the provider, calls
// `rules::decompose_delins`, and rebuilds N variants when the
// decomposition fires. For variants that don't decompose (most cases), the
// helper returns `vec![input]` and is effectively a no-op.

/// Per-coord-system-aware extraction of `(hgvs_start, hgvs_end, alt_bytes)`
/// from a variant whose edit is a literal `Delins` at simple positions
/// (no offsets, no uncertainty). Returns `None` for any variant shape that
/// can't be decomposed by the post-canonicalization split rules
/// (issues #160 / #165): non-Delins, intronic, uncertain boundary,
/// non-literal insert, etc.
/// Return the transcript-axis accession string for any `HgvsVariant`
/// kind that `apply_canonical_split` operates on. Used to build the
/// `CanonicalSplitSkipped` warning message — best-effort; returns
/// `<unknown>` for variant kinds that don't carry a single accession.
fn variant_accession_string(variant: &HgvsVariant) -> String {
    match variant {
        HgvsVariant::Genome(v) => v.accession.transcript_accession(),
        HgvsVariant::Cds(v) => v.accession.transcript_accession(),
        HgvsVariant::Tx(v) => v.accession.transcript_accession(),
        HgvsVariant::Rna(v) => v.accession.transcript_accession(),
        HgvsVariant::Mt(v) => v.accession.transcript_accession(),
        _ => "<unknown>".to_string(),
    }
}

fn extract_simple_delins(variant: &HgvsVariant) -> Option<(u64, u64, Vec<u8>)> {
    let (start, end, edit) = match variant {
        HgvsVariant::Genome(v) => simple_genome_loc_edit(&v.loc_edit)?,
        HgvsVariant::Cds(v) => simple_cds_loc_edit(&v.loc_edit)?,
        HgvsVariant::Tx(v) => simple_tx_loc_edit(&v.loc_edit)?,
        HgvsVariant::Rna(v) => simple_rna_loc_edit(&v.loc_edit)?,
        HgvsVariant::Mt(v) => simple_genome_loc_edit(&v.loc_edit)?,
        _ => return None,
    };
    let NaEdit::Delins { sequence, .. } = edit else {
        return None;
    };
    let InsertedSequence::Literal(seq) = sequence else {
        return None;
    };
    let alt: Vec<u8> = seq.bases().iter().map(|b| b.to_u8()).collect();
    Some((start, end, alt))
}

fn simple_genome_loc_edit(
    le: &LocEdit<Interval<GenomePos>, NaEdit>,
) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    let s = simple_genome_pos(le.location.start.as_single()?)?;
    let e = simple_genome_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_genome_pos(mu: &Mu<GenomePos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_special() || p.offset.is_some() {
        return None;
    }
    Some(p.base)
}

fn simple_cds_loc_edit(le: &LocEdit<Interval<CdsPos>, NaEdit>) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    // Only handle simple positive CDS positions (no UTR, no intronic, no
    // uncertainty). UTR delins decomposition would need its own coord-axis
    // logic and is out of scope for this fix.
    let s = simple_cds_pos(le.location.start.as_single()?)?;
    let e = simple_cds_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_cds_pos(mu: &Mu<CdsPos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_unknown() || p.is_intronic() || p.is_3utr() || p.base <= 0 {
        return None;
    }
    Some(p.base as u64)
}

fn simple_tx_loc_edit(le: &LocEdit<Interval<TxPos>, NaEdit>) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    let s = simple_tx_pos(le.location.start.as_single()?)?;
    let e = simple_tx_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_tx_pos(mu: &Mu<TxPos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_intronic() || p.is_downstream() || p.base <= 0 {
        return None;
    }
    Some(p.base as u64)
}

fn simple_rna_loc_edit(le: &LocEdit<Interval<RnaPos>, NaEdit>) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    let s = simple_rna_pos(le.location.start.as_single()?)?;
    let e = simple_rna_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_rna_pos(mu: &Mu<RnaPos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_intronic() || p.is_3utr() || p.base <= 0 {
        return None;
    }
    Some(p.base as u64)
}

/// Build a single HgvsVariant matching `template`'s coord-system kind /
/// accession / gene_symbol, with a new `[start_1based, end_1based]` location
/// and the given edit. Used by `build_split_variants` to spread the output
/// of `decompose_delins` back into a sequence of HgvsVariants.
fn build_variant_at(
    template: &HgvsVariant,
    start_1based: u64,
    end_1based: u64,
    edit: NaEdit,
) -> HgvsVariant {
    match template {
        HgvsVariant::Genome(g) => HgvsVariant::Genome(GenomeVariant {
            accession: g.accession.clone(),
            gene_symbol: g.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(start_1based), GenomePos::new(end_1based)),
                edit,
            ),
        }),
        HgvsVariant::Cds(c) => HgvsVariant::Cds(CdsVariant {
            accession: c.accession.clone(),
            gene_symbol: c.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(
                    CdsPos::new(start_1based as i64),
                    CdsPos::new(end_1based as i64),
                ),
                edit,
            ),
        }),
        HgvsVariant::Tx(t) => HgvsVariant::Tx(TxVariant {
            accession: t.accession.clone(),
            gene_symbol: t.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(
                    TxPos::new(start_1based as i64),
                    TxPos::new(end_1based as i64),
                ),
                edit,
            ),
        }),
        HgvsVariant::Rna(r) => HgvsVariant::Rna(RnaVariant {
            accession: r.accession.clone(),
            gene_symbol: r.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(
                    RnaPos::new(start_1based as i64),
                    RnaPos::new(end_1based as i64),
                ),
                edit,
            ),
        }),
        HgvsVariant::Mt(m) => HgvsVariant::Mt(MtVariant {
            accession: m.accession.clone(),
            gene_symbol: m.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(start_1based), GenomePos::new(end_1based)),
                edit,
            ),
        }),
        _ => unreachable!("build_variant_at called with non-NaEdit variant kind"),
    }
}

/// Build N HgvsVariants from a Vec<DelinsSubedit>. Position offsets in the
/// subedits are 0-indexed into the (per-variant-sized) ref slice; absolute
/// 1-based HGVS positions are recovered as `offset + hgvs_start`, where
/// `hgvs_start` is the variant's HGVS start position.
///
/// Spec rules implemented (see `general.md`, `substitution.md`):
///
/// 1. **Codon-frame exception** (`general.md:35-38`, issue #79 / #165).
///    When `codon_frame_aware` is true, the scan looks ahead at each
///    position for a `[Sub@i; Identity@i+1; Sub@i+2]` triplet whose CDS
///    endpoints (`hgvs_start + i`, `hgvs_start + i + 2`) share a codon.
///    Such a triplet emits as a single 3-base `delins` with alt
///    sequence `[Sub@i.alt, Identity@i+1.base, Sub@i+2.alt]`. The flag
///    is true only for `c.` (CDS) variants — `g.`, `n.`, `r.`, and `m.`
///    have no codon-frame and skip this branch. The exception is
///    deliberately narrow (length-3, exact pattern, in-codon endpoints)
///    so it matches the spec text "two variants separated by one
///    nucleotide, together affecting one amino acid".
///
/// 2. **Adjacent-substitution coalescence** (`substitution.md`,
///    issue #182). Consecutive `Substitution` sub-edits whose positions
///    are strictly adjacent (no gap, no `Inversion` or `IdentityAt`
///    between them) group into a single `delins` variant — "changes
///    involving two or more consecutive nucleotides are described as
///    deletion/insertion".
///
/// 3. **Inversion as a hard barrier** (issue #166). An `Inversion`
///    always emits standalone and breaks any in-flight substitution
///    run, preserving the inv-priority decomposition.
///
/// Singleton sub-runs stay as `Substitution`. `IdentityAt` not consumed
/// by a codon-frame triplet drops (an unchanged base is not an edit) and
/// always ends any in-flight substitution run — the gap means the
/// surrounding subs are no longer "consecutive".
fn build_split_variants(
    template: &HgvsVariant,
    subedits: Vec<DelinsSubedit>,
    hgvs_start: u64,
    codon_frame_aware: bool,
) -> Vec<HgvsVariant> {
    let abs = |idx: usize| -> u64 { idx as u64 + hgvs_start };

    let mut output: Vec<HgvsVariant> = Vec::new();
    // Pending run of strictly-adjacent Substitution sub-edits, in
    // left-to-right order. Each entry is `(position, reference, alternative)`
    // with `position` the 0-indexed offset into the variant's ref window.
    let mut run: Vec<(usize, Base, Base)> = Vec::new();

    let n = subedits.len();
    let mut i = 0;
    while i < n {
        // Codon-frame triplet lookahead: try to consume `[Sub; Identity; Sub]`
        // at offsets `[i, i+1, i+2]` whose endpoints share a codon. Only
        // fires for CDS variants (`codon_frame_aware`) and is the post-merge
        // half of issue #79: a pair of in-codon SNVs separated by one
        // unchanged base must render as a 3-base `delins`, even when the
        // pair sits inside a longer decomposition.
        if codon_frame_aware && i + 2 < n {
            if let (
                DelinsSubedit::Substitution {
                    position: p1,
                    alternative: a1,
                    ..
                },
                DelinsSubedit::IdentityAt {
                    position: pm,
                    base: bm,
                },
                DelinsSubedit::Substitution {
                    position: p3,
                    alternative: a3,
                    ..
                },
            ) = (&subedits[i], &subedits[i + 1], &subedits[i + 2])
            {
                if *pm == *p1 + 1 && *p3 == *p1 + 2 {
                    let cds_p1 = abs(*p1) as i64;
                    let cds_p3 = abs(*p3) as i64;
                    if merge::same_codon(cds_p1, cds_p3) {
                        // Codon-frame triplet preserved as a 3-base
                        // delins. `bm` is the unchanged ref byte from
                        // `decompose_delins`. Both `c.` and `r.` flow
                        // through this branch (issue #275 item 1 set
                        // `codon_frame_aware = true` for
                        // `HgvsVariant::Cds` and `HgvsVariant::Rna`),
                        // but no T/U recovery is needed here: the
                        // emitted delins is rendered by the per-variant
                        // formatter, and the `r.` formatter lowercases
                        // all of its bases (T → u included) when it
                        // prints the alt sequence. Forwarding the raw
                        // ref byte from `decompose_delins` is therefore
                        // safe for both coordinate systems.
                        flush_substitution_run(&mut output, template, hgvs_start, &mut run);
                        let s = abs(*p1);
                        let e = abs(*p3);
                        let alt_bases = vec![*a1, *bm, *a3];
                        output.push(build_variant_at(
                            template,
                            s,
                            e,
                            NaEdit::Delins {
                                sequence: InsertedSequence::Literal(Sequence::new(alt_bases)),
                                deleted: None,
                                deleted_length: None,
                            },
                        ));
                        i += 3;
                        continue;
                    }
                }
            }
        }

        match &subedits[i] {
            DelinsSubedit::Substitution {
                position,
                reference,
                alternative,
            } => {
                let breaks_run = matches!(run.last(), Some((prev, _, _)) if *prev + 1 != *position);
                if breaks_run {
                    flush_substitution_run(&mut output, template, hgvs_start, &mut run);
                }
                run.push((*position, *reference, *alternative));
            }
            DelinsSubedit::Inversion { start, end } => {
                flush_substitution_run(&mut output, template, hgvs_start, &mut run);
                // Half-open 0-indexed [start, end) of length L=end-start.
                // HGVS inclusive interval covers L bases starting at
                // abs(start) and ending at abs(start)+L-1 = abs(end-1).
                let s = abs(*start);
                let e = abs(*end) - 1;
                output.push(build_variant_at(
                    template,
                    s,
                    e,
                    NaEdit::Inversion {
                        sequence: None,
                        length: None,
                    },
                ));
            }
            // Drop IdentityAt: an unchanged base is not an edit. Outside of
            // the codon-frame triplet branch above, identities here are
            // either codon-frame-merge interior bases (issue #79) that did
            // not pair into a same-codon triplet, or outer bases absorbed
            // by `shorten_inversion`. An identity also ends any in-flight
            // substitution run — the gap means the surrounding subs are no
            // longer "consecutive".
            DelinsSubedit::IdentityAt { .. } => {
                flush_substitution_run(&mut output, template, hgvs_start, &mut run);
            }
        }
        i += 1;
    }
    flush_substitution_run(&mut output, template, hgvs_start, &mut run);
    output
}

/// Flush a pending run of consecutive substitution sub-edits into `output`.
/// A length-1 run emits a `Substitution`; a length-2+ run emits a single
/// `Delins` over `[run.first.position, run.last.position]` with `sequence`
/// = concatenated `alternative` bases. See `build_split_variants` for the
/// spec rationale (issue #182).
fn flush_substitution_run(
    output: &mut Vec<HgvsVariant>,
    template: &HgvsVariant,
    hgvs_start: u64,
    run: &mut Vec<(usize, Base, Base)>,
) {
    if run.is_empty() {
        return;
    }
    let abs = |idx: usize| -> u64 { idx as u64 + hgvs_start };
    if run.len() == 1 {
        let (position, reference, alternative) = run.pop().unwrap();
        let p = abs(position);
        output.push(build_variant_at(
            template,
            p,
            p,
            NaEdit::Substitution {
                reference,
                alternative,
            },
        ));
        return;
    }
    let s = abs(run.first().unwrap().0);
    let e = abs(run.last().unwrap().0);
    let alt_bases: Vec<Base> = run.drain(..).map(|(_, _, a)| a).collect();
    output.push(build_variant_at(
        template,
        s,
        e,
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(alt_bases)),
            deleted: None,
            deleted_length: None,
        },
    ));
}

/// Normalize DNA `T` and RNA `U` to a single byte (`T`) so byte-wise
/// comparison works across coord systems. Used by `apply_canonical_split`
/// to make every decomposition scan (rev-comp inv detection, per-position
/// sub / identity classification) T/U-agnostic for `r.` variants whose
/// alt bytes contain `U` while the transcript ref contains `T`.
fn normalize_t_u(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'U' => b'T',
            b'u' => b't',
            other => other,
        })
        .collect()
}

/// If `variants` has 1 element return it directly; if >1 wrap in a cis Allele.
fn wrap_allele_if_split(mut variants: Vec<HgvsVariant>) -> HgvsVariant {
    debug_assert!(!variants.is_empty(), "wrap_allele_if_split: empty input");
    if variants.len() == 1 {
        variants.pop().unwrap()
    } else {
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Cis))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    #[test]
    fn test_normalize_substitution_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();

        // Substitutions should not change
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_with_config() {
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime);
        let normalizer = Normalizer::with_config(provider, config);

        assert_eq!(
            normalizer.config().shuffle_direction,
            ShuffleDirection::FivePrime
        );
    }

    #[test]
    fn test_normalizer_handles_missing_transcript() {
        let provider = MockProvider::new(); // Empty provider
        let normalizer = Normalizer::new(provider);

        // Should return variant unchanged when transcript not found
        let variant = parse_hgvs("NM_MISSING.1:c.100del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
        // Verify output equals input (unchanged)
        assert_eq!(
            format!("{}", variant),
            format!("{}", result.unwrap()),
            "Missing transcript should return variant unchanged"
        );
    }

    #[test]
    fn test_normalize_deletion_shifts_3prime() {
        // NM_001234.1 has G repeat spanning exon boundaries
        // Exon 1: c.1-11, Exon 2: c.12-26, Exon 3: c.27+
        // G repeat is at c.9-c.33, but shift stops at exon boundary
        // c.10 is in exon 1 (ends at c.11), so deletion shifts to c.11
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.10del").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Should remain a deletion
        assert!(
            output.contains("del"),
            "Deletion should remain a deletion, got: {}",
            output
        );
        // Should shift from c.10 to c.11 (3'-most within exon 1)
        assert!(
            output.contains("c.11del"),
            "Deletion should shift 3' to exon boundary (c.11), got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_insertion_becomes_dup() {
        // NM_001234.1 has G repeat at CDS positions c.9-c.33
        // Inserting G after position 10 should shift 3' and become dup
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.10_11insG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Inserting G in a G-repeat should become dup
        assert!(
            output.contains("dup"),
            "Insertion of matching base should become dup, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_duplication_shifts_3prime() {
        // NM_001234.1 has G repeat spanning positions c.9-33 (25 G's)
        // Single-base duplications stay as simple dups (only 2+ base dups become repeat notation)
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.10dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' and stay as dup
        // c.10dup in GGGGG...GGG tract shifts to rightmost position (c.33) but stays as dup
        assert!(
            output.contains("dup"),
            "Single-base duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_delins_unchanged() {
        // A delins that doesn't simplify should stay as delins
        // Deleting G and inserting AT is not a dup pattern
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delinsAT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("delinsAT"),
            "Delins should remain unchanged, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_single_base_delins_becomes_substitution() {
        // HGVS edit-type priority: a 1→1 delins with ref!=alt must be expressed
        // as a substitution. Transcript NM_000088.3 starts ATGCCCAAGG...; position
        // 10 is G. c.10delinsT replaces G with T → c.10G>T.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delinsT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10G>T");
    }

    #[test]
    fn test_normalize_single_base_delins_same_base_becomes_identity() {
        // Per HGVS, a delins whose insert equals the reference is identity (=).
        // Transcript NM_000088.3 starts ATGCCCAAGG...; position 10 is G.
        // c.10delinsG produces no change → c.10=.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delinsG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10=");
    }

    #[test]
    fn test_normalize_multi_base_delins_same_seq_becomes_identity() {
        // Transcript NM_000088.3 starts ATG at positions 1-3.
        // c.1_3delinsATG produces no change → c.1_3=.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.1_3delinsATG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.1_3=");
    }

    #[test]
    fn test_normalize_multi_base_delete_delins_to_pure_deletion() {
        // c.10_11delinsT against NM_000088.3 (c.10_11 = GT). The shared `T`
        // suffix consumes the inserted base entirely, leaving a single-base
        // deletion at c.10. Per HGVS minimal-form rules (sub > del > inv >
        // dup > ins) the canonical output is a pure deletion, not a delins.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10_11delinsT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10del");
    }

    #[test]
    fn test_normalize_empty_insert_delins_becomes_deletion() {
        // HGVS spec: a delins whose inserted sequence is empty is semantically
        // a deletion and must be rendered as `del`. Issue #81 item A3.
        // Transcript NM_000088.3 starts ATGCCCAAGG…; c.10delins (empty insert)
        // is a deletion of position 10 → c.10del.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delins").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10del");
    }

    #[test]
    fn test_normalize_empty_insert_multi_base_delins_becomes_deletion() {
        // Multi-base form: c.10_11delins (deletes "GT" at positions 10-11,
        // inserts nothing) → del, then the spec's 3'-rule shifts the deletion
        // to c.11_12del because ref[10]=G == ref[12]=G. Issue #81 item A3.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10_11delins").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.11_12del");
    }

    #[test]
    fn test_normalize_delins_to_dup_still_works() {
        // Regression guard: adding identity/substitution checks before the dup
        // check must not block legitimate dup conversions. ref[5] = C;
        // c.5delinsCC matches the dup pattern (insert == deleted twice) and
        // must still normalize to dup.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.5delinsCC").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);
        assert!(
            output.contains("dup"),
            "delins matching dup pattern should normalize to dup, got: {}",
            output
        );
        assert!(
            !output.contains("delins"),
            "delins should not survive the dup conversion, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_delins_different_bases_becomes_substitution() {
        // c.1_3delinsACG against NM_000088.3 (c.1_3 = ATG). The shared `A`
        // prefix and `G` suffix collapse the delins to T -> C at c.2 per the
        // HGVS minimal-form rule (sub > del > inv > dup > ins).
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.1_3delinsACG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.2T>C");

        // r. twin: same collapse on the RNA path. NM_000088.3 has cds_start = 1,
        // so r.1_3 maps to the same ATG run and the residual is at r.2. The
        // residual ref byte is DNA `T` (sourced from the c. axis), which the
        // RNA `Display` path canonicalizes to `u` per the HGVS RNA alphabet
        // (see issue #276 / `NaEdit::Substitution` in src/hgvs/edit.rs).
        let variant = parse_hgvs("NM_000088.3:r.1_3delinsacg").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:r.2u>c");
    }

    #[test]
    fn test_normalize_substitution_ref_equals_alt_becomes_identity() {
        // Per HGVS, a substitution where the reference and alternative bases
        // are identical produces no change and must be expressed using
        // identity notation (`=`). SNV companion to the same-base delins rule.
        // Transcript NM_000088.3 starts ATGCCCAAGG...; position 10 is G.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10G>G").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10=");
    }

    #[test]
    fn test_normalize_substitution_ref_equals_alt_first_position() {
        // Boundary check: the rule must fire at position 1 (first base).
        // Transcript NM_000088.3 starts ATG...; position 1 is A.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.1A>A").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.1=");
    }

    #[test]
    fn test_normalize_substitution_ref_not_equal_alt_unchanged() {
        // Regression guard: a real SNV must not be rewritten to identity.
        // c.10G>T at position 10 (G) is a valid substitution.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10G>T").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10G>T");
    }

    #[test]
    fn test_normalize_substitution_ref_equals_alt_without_provider_data() {
        // The A4 rule is purely syntactic, so it must fire even when the
        // provider has no transcript loaded — matching the spec's stance that
        // `c.123C>C` is "not allowed" regardless of reference availability.
        // Spec example: docs/recommendations/DNA/other.md (HGVS v21.0).
        let provider = MockProvider::new(); // empty — no transcripts
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_004006.2:c.123C>C").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_004006.2:c.123=");
    }

    #[test]
    fn test_normalize_protein_substitution_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Protein substitution variants should pass through unchanged
        let variant = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_protein_deletion_removes_redundant_sequence() {
        // Redundant sequence removal: p.Val600delVal → p.Val600del
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600delVal").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        // Should remove redundant "Val" from the deletion
        assert_eq!(
            output, "NP_000079.2:p.Val600del",
            "Redundant sequence should be removed from deletion"
        );
    }

    #[test]
    fn test_normalize_protein_deletion_range_removes_redundant_sequence() {
        // Redundant sequence removal for range: p.Lys23_Glu25delLysAlaGlu → p.Lys23_Glu25del
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Lys23_Glu25delLysAlaGlu").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        // Should remove redundant sequence from the deletion
        assert_eq!(
            output, "NP_000079.2:p.Lys23_Glu25del",
            "Redundant sequence should be removed from range deletion"
        );
    }

    #[test]
    fn test_normalize_protein_deletion_non_matching_sequence_unchanged() {
        // Non-matching sequence should stay: p.Val600delGlu should NOT change
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600delGlu").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        // Should NOT remove non-matching sequence
        assert_eq!(
            output, "NP_000079.2:p.Val600delGlu",
            "Non-matching sequence should not be removed"
        );
    }

    #[test]
    fn test_normalize_protein_deletion_no_sequence_unchanged() {
        // Deletion without sequence should stay unchanged
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600del").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        assert_eq!(
            output, "NP_000079.2:p.Val600del",
            "Deletion without sequence should remain unchanged"
        );
    }

    #[test]
    fn test_normalize_protein_duplication_unchanged() {
        // Duplications should pass through unchanged
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600dup").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_protein_frameshift_unchanged() {
        // Frameshifts should pass through unchanged
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Arg97ProfsTer23").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_protein_reference_validation_match() {
        // Test that validation passes when amino acid matches reference
        // NP_TEST.1 has: M at position 1, V at position 2
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position 1 = M (Met), Position 2 = V (Val)
        let variant = parse_hgvs("NP_TEST.1:p.Met1Val").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Validation should pass for matching amino acid"
        );
    }

    #[test]
    fn test_normalize_protein_reference_validation_mismatch() {
        // Test that validation fails when amino acid doesn't match reference
        // NP_TEST.1 has M at position 1, but we claim it's Val
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position 1 is M (Met), not V (Val)
        let variant = parse_hgvs("NP_TEST.1:p.Val1Glu").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_err(),
            "Validation should fail for mismatched amino acid"
        );

        if let Err(crate::error::FerroError::AminoAcidMismatch {
            position,
            expected,
            found,
            ..
        }) = result
        {
            assert_eq!(position, 1);
            assert_eq!(expected, "Val");
            assert_eq!(found, "M");
        } else {
            panic!("Expected AminoAcidMismatch error");
        }
    }

    #[test]
    fn test_normalize_protein_reference_validation_deletion() {
        // Test validation for deletion variants
        // NP_TEST.1 has V at position 2
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position 2 = V (Val) - should pass
        let variant = parse_hgvs("NP_TEST.1:p.Val2del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Validation should pass for matching deletion position"
        );
    }

    #[test]
    fn test_normalize_protein_reference_validation_missing_protein() {
        // Test that missing protein data skips validation (doesn't fail)
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // NP_MISSING.1 doesn't exist in provider
        let variant = parse_hgvs("NP_MISSING.1:p.Val600Glu").unwrap();
        let result = normalizer.normalize(&variant);
        // Should NOT error - just skip validation when protein not available
        assert!(
            result.is_ok(),
            "Missing protein should skip validation, not fail"
        );
    }

    #[test]
    fn test_normalize_rna_substitution_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // RNA substitutions should pass through unchanged
        let variant = parse_hgvs("NM_000088.3:r.10a>g").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_rna_deletion_shifts_3prime() {
        // NM_001234.1 has G repeat at positions 13-37 (position 37 is actually T,
        // but the normalization shifts to the 3'-most position)
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.14del").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Should remain a deletion and shift 3'
        assert!(
            output.contains("del"),
            "Deletion should remain a deletion, got: {}",
            output
        );
        // Verify it shifted 3' (position should be > 14)
        assert!(
            output.contains("r.37del"),
            "Deletion should shift 3' to position 37, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_rna_insertion_becomes_dup() {
        // NM_001234.1 has G repeat at positions 13-36
        // Inserting g after position 14 should shift 3' and become dup
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.14_15insg").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Inserting g in a G-repeat should become dup
        assert!(
            output.contains("dup"),
            "Insertion of matching base should become dup, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_rna_duplication_shifts_3prime() {
        // NM_001234.1 has G repeat - single-base duplications stay as simple dups
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.14dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' and stay as dup
        assert!(
            output.contains("dup"),
            "Single-base RNA duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_rna_intronic_returns_error() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Intronic RNA variants should return an error
        let variant = parse_hgvs("NM_001234.1:r.10+5del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_err(), "Intronic RNA variant should return error");
    }

    #[test]
    fn test_normalize_rna_missing_transcript_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Missing transcript should return variant unchanged
        let variant = parse_hgvs("NM_MISSING.1:r.100del").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_null_allele() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Null alleles should pass through
        let variant = HgvsVariant::NullAllele;
        let normalized = normalizer.normalize(&variant).unwrap();
        assert!(matches!(normalized, HgvsVariant::NullAllele));
    }

    #[test]
    fn test_normalize_unknown_allele() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Unknown alleles should pass through
        let variant = HgvsVariant::UnknownAllele;
        let normalized = normalizer.normalize(&variant).unwrap();
        assert!(matches!(normalized, HgvsVariant::UnknownAllele));
    }

    #[test]
    fn test_normalize_allele_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Allele variants should normalize each component
        // Using substitutions which should remain unchanged
        let variant = parse_hgvs("[NM_000088.3:c.10A>G;NM_000088.3:c.20C>T]").unwrap();
        let result = normalizer.normalize(&variant).unwrap();

        // Verify it's still an allele
        assert!(matches!(result, HgvsVariant::Allele(_)));

        // Verify output is the canonical compact form (ACC:c.[edit1;edit2])
        assert_eq!(
            format!("{}", result),
            "NM_000088.3:c.[10A>G;20C>T]",
            "Allele display should use canonical compact form"
        );
    }

    #[test]
    fn test_normalize_5prime_direction() {
        // Test that 5' direction shifts toward 5' end instead of 3'
        // NM_001234.1 has G repeat spanning exons
        // Exon 2: c.12-c.26
        // With 5' direction, deletion at c.20 should shift toward c.12
        // Note: Actual shift depends on boundary handling
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime);
        let normalizer = Normalizer::with_config(provider, config);

        // Delete G at position 20 (middle of G repeat in exon 2)
        let variant = parse_hgvs("NM_001234.1:c.20del").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("del"),
            "Should remain a deletion, got: {}",
            output
        );
        // With 5' direction within exon 2, should shift toward 5' boundary
        // The exact position depends on exon boundary handling
        assert!(
            output.contains("c.13del") || output.contains("c.12del"),
            "5' direction should shift deletion toward exon start, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_3prime_direction() {
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime);
        let normalizer = Normalizer::with_config(provider, config);

        let variant = parse_hgvs("NM_000088.3:c.10del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_with_cross_boundaries() {
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().allow_crossing_boundaries();
        let normalizer = Normalizer::with_config(provider, config);

        let variant = parse_hgvs("NM_000088.3:c.10del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_genomic_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Genomic variants with missing sequence should pass through unchanged
        let variant = parse_hgvs("NC_000001.11:g.12345del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_tx_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NR_000001.1:n.100del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_mt_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // MT variants pass through unchanged
        let variant = parse_hgvs("NC_012920.1:m.100A>G").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert!(matches!(normalized, HgvsVariant::Mt(_)));
    }

    #[test]
    fn test_cds_to_tx_pos_utr5() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // 5' UTR position (negative)
        let cds_pos = CdsPos {
            base: -5,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&cds_pos, 10, Some(50));
        assert!(result.is_ok());
    }

    #[test]
    fn test_cds_to_tx_pos_utr3() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // 3' UTR position
        let cds_pos = CdsPos {
            base: 5,
            offset: None,
            utr3: true,
        };
        let result = normalizer.cds_to_tx_pos(&cds_pos, 10, Some(50));
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 55);
    }

    #[test]
    fn test_cds_to_tx_pos_coding() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Normal coding position
        let cds_pos = CdsPos {
            base: 10,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&cds_pos, 5, Some(50));
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 14); // 5 + 10 - 1 = 14
    }

    #[test]
    fn test_tx_to_cds_pos_utr5() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position before CDS start
        let result = normalizer.tx_to_cds_pos(3, 10, Some(50));
        assert!(result.is_ok());
        let cds_pos = result.unwrap();
        assert!(cds_pos.base < 0);
        assert!(!cds_pos.utr3);
    }

    #[test]
    fn test_tx_to_cds_pos_utr3() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position after CDS end
        let result = normalizer.tx_to_cds_pos(55, 10, Some(50));
        assert!(result.is_ok());
        let cds_pos = result.unwrap();
        assert!(cds_pos.utr3);
    }

    #[test]
    fn test_tx_to_cds_pos_coding() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Normal coding position
        let result = normalizer.tx_to_cds_pos(20, 10, Some(50));
        assert!(result.is_ok());
        let cds_pos = result.unwrap();
        assert!(!cds_pos.utr3);
        assert_eq!(cds_pos.base, 11); // 20 - 10 + 1 = 11
    }

    #[test]
    fn test_config_default() {
        let config = NormalizeConfig::default();
        assert_eq!(config.shuffle_direction, ShuffleDirection::ThreePrime);
        assert!(!config.cross_boundaries);
    }

    #[test]
    #[allow(deprecated)]
    fn test_config_builder() {
        let config = NormalizeConfig::default()
            .with_direction(ShuffleDirection::FivePrime)
            .allow_crossing_boundaries()
            .skip_validation();

        assert_eq!(config.shuffle_direction, ShuffleDirection::FivePrime);
        assert!(config.cross_boundaries);
        // skip_validation now sets RefSeqMismatch to SilentCorrect
        assert!(!config.should_reject_ref_mismatch());
        assert!(!config.should_warn_ref_mismatch());
    }

    #[test]
    fn test_duplication_3prime_shift_two_bases() {
        // Test the exact scenario from ClinVar: c.4159dup vs c.4160dup
        // When duplicating an A in a homopolymer tract (AA),
        // HGVS requires converting to repeat notation.
        //
        // NM_888888.1 sequence: ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        // Positions (1-based):  12345678901234567890...
        // c.8 = A, c.9 = A (the "AA" in "GAA")
        //
        // Single-base duplications stay as simple dups
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_888888.1:c.8dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' and stay as dup
        assert!(
            output.contains("dup"),
            "Single-base duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_duplication_3prime_shift_three_bases() {
        // Test with three consecutive identical bases (TTT)
        //
        // NM_888888.1 sequence: ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        // Positions:           ...              2021222324...
        // c.20 = G, c.21 = T, c.22 = T, c.23 = T, c.24 = G
        //
        // Single-base duplications stay as simple dups
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_888888.1:c.21dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' to c.23 and stay as dup
        assert!(
            output.contains("dup"),
            "Single-base duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_duplication_no_shift_when_unique() {
        // Test that a duplication of a unique base doesn't shift
        //
        // NM_888888.1 sequence: ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        // c.7 = G (followed by AA, so no G to shift to)
        //
        // c.7dup should stay as c.7dup
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_888888.1:c.7dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("dup"),
            "Should remain a duplication, got: {}",
            output
        );
        assert!(
            output.contains("c.7dup"),
            "Duplication of unique G at c.7 should not shift, got: {}",
            output
        );
    }

    // =====================================================================
    // Exon-Intron Boundary Spanning Variant Tests
    // =====================================================================

    /// Create a provider with a transcript that has genomic coordinates and introns
    /// for testing boundary-spanning variant normalization.
    ///
    /// Transcript structure (NM_BOUNDARY.1):
    /// - Gene: BOUNDARY
    /// - Strand: Plus
    /// - Chromosome: chr1
    ///
    /// Genomic layout (chr1):
    /// Position: 1000                   1020  1030                   1050  1060                   1080
    ///           |-------- Exon 1 ------|      |-------- Exon 2 ------|      |-------- Exon 3 ------|
    ///           ATGCCCAAAGGGTTTAGGCCC       AAAGGGTTTAGGCCCAAAAAA       GGGTTTAGGCCCAAATGA
    ///                                 ^^^  ^^^                   ^^^  ^^^
    ///                              intron 1                   intron 2
    ///
    /// Transcript positions (tx):
    /// - Exon 1: tx 1-20 = genomic 1000-1019
    /// - Intron 1: genomic 1020-1029 (10 bp)
    /// - Exon 2: tx 21-40 = genomic 1030-1049
    /// - Intron 2: genomic 1050-1059 (10 bp)
    /// - Exon 3: tx 41-58 = genomic 1060-1077
    ///
    /// CDS: starts at tx 1 (no 5' UTR for simplicity)
    /// CDS positions: c.1 = tx 1, c.20 = tx 20 (last of exon 1), c.21 = tx 21 (first of exon 2)
    ///
    /// Intron 1 sequence (g.1020-1029): "GTAAGCTAGG" (10 bp)
    ///   - c.20+1 = g.1020 (G)
    ///   - c.20+10 = g.1029 (G)
    ///   - c.21-10 = g.1020 (G)
    ///   - c.21-1 = g.1029 (G)
    ///
    /// Intron 2 sequence (g.1050-1059): "GTAAGTAAGG" (10 bp)
    fn make_boundary_test_provider() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Build transcript sequence (exons only, spliced)
        // Exon 1 (20bp): ATGCCCAAAGGGTTTAGGCC (ends with CC at exon boundary)
        // Exon 2 (20bp): AAAGGGTTTAGGCCCAAAAA (AA repeat at both boundaries)
        // Exon 3 (18bp): GGGTTTAGGCCCAAATGA
        // Total: 58bp transcript
        let tx_seq = "ATGCCCAAAGGGTTTAGGCCAAAGGGTTTAGGCCCAAAAAGGGTTTAGGCCCAAATGA";

        // Build genomic sequence around the transcript
        // We'll create 100bp before, the gene region, and 100bp after
        // Gene region: exon1 + intron1 + exon2 + intron2 + exon3
        // = 20 + 10 + 20 + 10 + 18 = 78bp
        //
        // Genomic: 900-999 (padding) + 1000-1019 (exon1) + 1020-1029 (intron1) +
        //          1030-1049 (exon2) + 1050-1059 (intron2) + 1060-1077 (exon3) + 1078+ (padding)
        let mut genomic_seq = String::new();

        // Padding before (positions 0-999, 1000 bytes at 0-based)
        for _ in 0..1000 {
            genomic_seq.push('N');
        }

        // Exon 1 (positions 1000-1019, 0-based 1000-1019)
        genomic_seq.push_str("ATGCCCAAAGGGTTTAGGCC");

        // Intron 1 (positions 1020-1029) - with splice consensus
        // Note: The intron has AAA at the end (1027-1029) to test shifting
        genomic_seq.push_str("GTAAGCTAAA");

        // Exon 2 (positions 1030-1049)
        genomic_seq.push_str("AAAGGGTTTAGGCCCAAAAA");

        // Intron 2 (positions 1050-1059) - with AAA at start for testing
        genomic_seq.push_str("AAAGTAAGGG");

        // Exon 3 (positions 1060-1077)
        genomic_seq.push_str("GGGTTTAGGCCCAAATGA");

        // Padding after (100 bytes)
        for _ in 0..100 {
            genomic_seq.push('N');
        }

        // Add genomic sequence to provider
        provider.add_genomic_sequence("chr1", genomic_seq);

        // Create transcript with exons that have genomic coordinates
        provider.add_transcript(Transcript {
            id: "NM_BOUNDARY.1".to_string(),
            gene_symbol: Some("BOUNDARY".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(58),
            exons: vec![
                Exon::with_genomic(1, 1, 20, 1000, 1019),
                Exon::with_genomic(2, 21, 40, 1030, 1049),
                Exon::with_genomic(3, 41, 58, 1060, 1077),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1077),
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        provider
    }

    /// Minus-strand mirror of `make_boundary_test_provider`.
    ///
    /// Same transcript sequence as the plus-strand fixture, so c.40_40+3
    /// still spans the same poly-A region in transcript view. The genomic
    /// content at the gene region is the reverse complement of each exon
    /// (so RC of the genomic plus strand recovers `tx_seq`), and the
    /// exon-to-genomic mapping is reversed: tx 1 maps to the high genomic
    /// end (g.1077) and tx 58 to the low end (g.1000). Intron 2 is laid
    /// out so that c.40+1..c.40+4 read as `A` in transcript view, putting
    /// the boundary-spanning dup inside the same 5-A tract that the plus
    /// fixture exercises.
    fn make_boundary_test_provider_minus() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        let tx_seq = "ATGCCCAAAGGGTTTAGGCCAAAGGGTTTAGGCCCAAAAAGGGTTTAGGCCCAAATGA";

        let mut genomic_seq = String::new();
        for _ in 0..1000 {
            genomic_seq.push('N');
        }
        // Exon 3 region (g.1000-1017): RC of tx[41..58] ("GGGTTTAGGCCCAAATGA").
        genomic_seq.push_str("TCATTTGGGCCTAAACCC");
        // Intron 2 (g.1018-1027): the last four bases (g.1024-1027) are 'T',
        // so c.40+1..c.40+4 read as 'A' in transcript view, extending the
        // exonic poly-A across the boundary.
        genomic_seq.push_str("AAAGTATTTT");
        // Exon 2 region (g.1028-1047): RC of tx[21..40] ("AAAGGGTTTAGGCCCAAAAA").
        genomic_seq.push_str("TTTTTGGGCCTAAACCCTTT");
        // Intron 1 (g.1048-1057): mirrors the plus fixture's intron 1 content.
        genomic_seq.push_str("GTAAGCTAAA");
        // Exon 1 region (g.1058-1077): RC of tx[1..20] ("ATGCCCAAAGGGTTTAGGCC").
        genomic_seq.push_str("GGCCTAAACCCTTTGGGCAT");
        for _ in 0..100 {
            genomic_seq.push('N');
        }

        provider.add_genomic_sequence("chr1", genomic_seq);

        provider.add_transcript(Transcript {
            id: "NM_BOUNDARYM.1".to_string(),
            gene_symbol: Some("BOUNDARY_M".to_string()),
            strand: Strand::Minus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(58),
            exons: vec![
                Exon::with_genomic(1, 1, 20, 1058, 1077),
                Exon::with_genomic(2, 21, 40, 1028, 1047),
                Exon::with_genomic(3, 41, 58, 1000, 1017),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1077),
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        provider
    }

    #[test]
    fn test_boundary_spanning_exonic_to_intronic_del() {
        // Test: c.20_20+3del - deletion from last exon base into intron
        // c.20 = last base of exon 1 (C at g.1019)
        // c.20+3 = 3rd intronic base (A at g.1022)
        // Deletes: C (exonic) + GTA (intronic) = 4 bases
        //
        // This should normalize without error (using genomic space)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+3del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("del"),
            "Should remain a deletion, got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_intronic_to_exonic_del() {
        // Test: c.21-3_23del - deletion from intron into exon
        // c.21-3 = 3rd base before exon 2 (A at g.1027)
        // c.23 = 3rd base of exon 2 (A at g.1032)
        // Deletes: AAA (intronic) + AAA (exonic) = 6 bases
        //
        // The intron ends with AAA (g.1027-1029) and exon starts with AAA (g.1030-1032)
        // This is a repeat, so normalization might shift
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.21-3_23del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("del"),
            "Should remain a deletion, got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_same_base_position() {
        // Test: c.20_20+5del - deletion starting and ending at same CDS base
        // Start is exonic (c.20), end is intronic (c.20+5)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+5del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Same-base boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_boundary_splice_site_plus1() {
        // Test: c.20_20+1del - deletion of last exon base + splice donor (GT)
        // This is a clinically important splice site variant
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+1del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Splice site +1 deletion should normalize, got error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_boundary_splice_site_minus1() {
        // Test: c.21-1_22del - deletion of splice acceptor + first exon bases
        // c.21-1 = last intronic base before exon 2 (A at g.1029)
        // c.22 = 2nd base of exon 2 (A at g.1031)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.21-1_22del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Splice site -1 deletion should normalize, got error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_boundary_spanning_dup() {
        // Test: c.40_40+3dup — duplication spanning exon-intron boundary,
        // landing on a 4-A poly-A region. Per HGVS spec (repeated.md):
        //
        //   > This restriction only applies to the coding sequence,
        //   > which does not include the introns or the UTR sequence.
        //
        // The codon-frame `unit_len % 3 == 0` restriction does NOT apply
        // to boundary-spanning variants (mixed exon/intron context is
        // not "purely coding sequence"), so `normalize_boundary_spanning_cds`
        // passes `is_coding=false` to `normalize_na_edit`. The dup must
        // therefore canonicalize to `A[N]` repeat notation, not the gated
        // `ins<literal>` fallback. (B4-remaining, issue #209.)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.40_40+3dup").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning duplication should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("A[") && !output.contains("ins"),
            "Boundary-spanning multi-copy dup spanning into intron must emit \
             `A[N]` repeat notation (intron exempts the codon-frame gate per \
             repeated.md), got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_dup_minus_strand() {
        // Minus-strand mirror of `test_boundary_spanning_dup`. Pins the
        // strand-specific flip in `normalize_boundary_spanning_cds`: the
        // genomic-strand window is RC of the transcript view, so without
        // flipping, repeat detection would inspect the wrong alphabet and
        // miss the A homopolymer. With `is_coding=false` (the
        // boundary-spanning context's intronic exemption — B4-remaining),
        // the canonical form is `A[N]` repeat notation on the transcript-
        // view bytes.
        let provider = make_boundary_test_provider_minus();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARYM.1:c.40_40+3dup").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning duplication should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("A[") && !output.contains("ins"),
            "Minus-strand boundary-spanning multi-copy dup must emit `A[N]` \
             repeat notation (intron exempts the codon-frame gate per \
             repeated.md), got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_del_emits_repeat_notation() {
        // Test: `c.40_40+3del` — deletion spanning exon-intron boundary
        // into the same 4-A poly-A region used by
        // `test_boundary_spanning_dup`. Symmetric counterpart on the
        // `del` branch (`deletion_to_repeat`) of the codon-frame gate
        // fix. Per the intronic exemption (issue #209 B4-remaining),
        // boundary-spanning context passes `is_coding=false`, so the
        // del of repeat-unit bases must canonicalize to `A[N-k]` over
        // the reference-tract extent rather than falling back to plain
        // `del`.
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.40_40+3del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("A[") && !output.contains("del"),
            "Boundary-spanning multi-copy del spanning into intron must \
             emit `A[N]` repeat notation (intron exempts the codon-frame \
             gate per repeated.md), got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_delins() {
        // Test: c.20_20+2delinsTTT - delins spanning exon-intron boundary
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+2delinsTTT").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning delins should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("delins") || output.contains(">"),
            "Should remain a delins or become substitution, got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_no_genomic_data_returns_error() {
        // Test that without genomic data, we still get the ExonIntronBoundary error
        let provider = MockProvider::with_test_data(); // No genomic data
        let normalizer = Normalizer::new(provider);

        // NM_001234.1 doesn't have genomic coordinates
        let variant = parse_hgvs("NM_001234.1:c.11_11+3del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_err(),
            "Boundary-spanning without genomic data should return error"
        );
    }

    #[test]
    fn test_deletion_3prime_shift_consecutive_bases() {
        // Test case simulating NM_001408491.1:c.517delA -> should become c.518del
        // Create a transcript with consecutive A's at positions that should shift
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Create a sequence where c.517 and c.518 are both 'A'
        // CDS starts at position 1, so c.N = transcript position N
        // Put "AA" at positions 517-518 (1-based)
        // Sequence: 516 bases of padding + "AA" + more padding
        let mut seq = String::new();
        for _ in 0..516 {
            seq.push('G'); // Padding (not A to ensure we see the shift)
        }
        seq.push('A'); // Position 517 (c.517)
        seq.push('A'); // Position 518 (c.518)
        for _ in 519..=600 {
            seq.push('G'); // More padding
        }

        provider.add_transcript(crate::reference::transcript::Transcript {
            id: "NM_777777.1".to_string(),
            gene_symbol: Some("SHIFTTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(seq.clone()),
            cds_start: Some(1),
            cds_end: Some(600),
            exons: vec![Exon::new(1, 1, 600)], // Single exon covering all
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        let normalizer = Normalizer::new(provider);

        // Parse c.517delA - deleting the first A
        let variant = parse_hgvs("NM_777777.1:c.517delA").unwrap();

        // Debug: print the sequence around positions 517-518
        println!("Sequence length: {}", seq.len());
        println!(
            "Position 516 (0-based 515): {}",
            seq.chars().nth(515).unwrap()
        );
        println!(
            "Position 517 (0-based 516): {}",
            seq.chars().nth(516).unwrap()
        );
        println!(
            "Position 518 (0-based 517): {}",
            seq.chars().nth(517).unwrap()
        );
        println!(
            "Position 519 (0-based 518): {}",
            seq.chars().nth(518).unwrap()
        );

        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        println!("Input:  NM_777777.1:c.517delA");
        println!("Output: {}", output);

        // The deletion should shift from c.517 to c.518 (3' rule)
        // because both positions are 'A'
        assert!(
            output.contains("c.518del"),
            "Deletion at c.517 should shift to c.518 (3' rule), got: {}",
            output
        );
    }

    #[test]
    fn test_deletion_3prime_shift_with_utr() {
        // Same test but with a 5' UTR (cds_start > 1)
        // This simulates real transcripts more accurately
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Create a transcript with 100bp 5' UTR
        // CDS starts at position 101, so:
        // c.1 = tx position 101
        // c.517 = tx position 617
        // c.518 = tx position 618
        let utr_len = 100;
        let mut seq = String::new();

        // 5' UTR (100 bases)
        for _ in 0..utr_len {
            seq.push('T');
        }
        // CDS: 516 bases of G padding, then "AA", then more G
        for _ in 0..516 {
            seq.push('G');
        }
        seq.push('A'); // tx position 617 = c.517
        seq.push('A'); // tx position 618 = c.518
        for _ in 0..100 {
            seq.push('G');
        }

        let seq_len = seq.len();
        println!("Test with UTR:");
        println!("  Sequence length: {}", seq_len);
        println!("  CDS start (1-based): 101");
        println!("  c.517 = tx position 617 (0-based 616)");
        println!("  c.518 = tx position 618 (0-based 617)");
        println!(
            "  tx pos 617 (0-based 616): {}",
            seq.chars().nth(616).unwrap()
        );
        println!(
            "  tx pos 618 (0-based 617): {}",
            seq.chars().nth(617).unwrap()
        );

        provider.add_transcript(crate::reference::transcript::Transcript {
            id: "NM_666666.1".to_string(),
            gene_symbol: Some("UTRTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(seq.clone()),
            cds_start: Some(101),
            cds_end: Some(seq_len as u64),
            exons: vec![Exon::new(1, 1, seq_len as u64)], // Single exon
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_666666.1:c.517delA").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        println!("Input:  NM_666666.1:c.517delA");
        println!("Output: {}", output);

        assert!(
            output.contains("c.518del"),
            "Deletion at c.517 should shift to c.518 (3' rule) even with UTR, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_inverted_range_insertion_no_panic() {
        // Regression: ClinVar pattern NC_000011.10:g.5238138_5153222insTATTT
        // has start > end (inverted range).  Previously caused a panic in
        // insertion_is_duplication due to slice index out of bounds.
        // The normalizer should return an error, not panic.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NC_000011.10:g.5238138_5153222insTATTT").unwrap();
        let result = normalizer.normalize(&variant);
        // It's fine if this returns Ok (unchanged) or Err (validation failure),
        // but it must NOT panic.
        let _ = result;
    }

    #[test]
    fn test_delins_should_not_shift() {
        // HGVS spec: delins should NOT be 3' shifted like del/dup/ins
        // This test ensures we don't incorrectly shift delins positions
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Create a transcript where delins could theoretically shift if we were wrong
        // Sequence: ...GGAATTCC... where we do c.10_11delinsXX
        // If incorrectly shifted, it might become c.11_12delinsXX
        let seq = "GGGGGGGGGGAATTCCGGGGGGGGGG".to_string(); // c.10=A, c.11=A, c.12=T, c.13=T

        provider.add_transcript(crate::reference::transcript::Transcript {
            id: "NM_555555.1".to_string(),
            gene_symbol: Some("DELINSTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(seq),
            cds_start: Some(1),
            cds_end: Some(26),
            exons: vec![Exon::new(1, 1, 26)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        let normalizer = Normalizer::new(provider);

        // Test delins - should NOT shift
        let variant = parse_hgvs("NM_555555.1:c.10_11delinsTT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("c.10_11delins"),
            "Delins should NOT be shifted (HGVS spec), got: {}",
            output
        );
    }

    #[test]
    fn test_cds_to_tx_pos_utr5_underflow() {
        // cds_start=5, base=-6 → 5 + (-6) - 1 = -2, should return Err not wrap to u64::MAX
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let pos = CdsPos {
            base: -6,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&pos, 5, Some(38));
        assert!(
            result.is_err(),
            "cds_to_tx_pos should return Err for positions before transcript start, got: {:?}",
            result
        );
    }

    #[test]
    fn test_cds_to_tx_pos_utr5_valid() {
        // HGVS numbering skips c.0, so c.-N maps to tx position
        // cds_start - N. For cds_start=5, c.-3 → tx = 5 + (-3) = 2.
        // Issue #97 — the previous formula `cds_start + base - 1`
        // double-counted the gap and returned tx 1 (the c.-4 base).
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let pos = CdsPos {
            base: -3,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&pos, 5, Some(38));
        assert_eq!(result.unwrap(), 2);
    }

    #[test]
    fn test_normalize_cds_utr5_deep_negative() {
        // A deeply negative 5' UTR position that would overflow should return an error, not panic
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        // c.-88 with cds_start=5 → 5 + (-88) - 1 = -84, which would wrap to huge u64
        let variant = parse_hgvs("NM_001234.1:c.-88A>G").unwrap();
        let result = normalizer.normalize(&variant);
        // The primary check is that this doesn't panic.
        let _ = result;
    }

    #[test]
    fn test_normalize_unknown_offset_returns_unchanged() {
        // Variants with ? offsets (sentinel values i64::MAX/MIN) should return unchanged
        // because we can't normalize with indeterminate boundaries
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // c.-85-?_834+?del has unknown offsets on both positions
        let variant = parse_hgvs("NM_000088.3:c.-85-?_834+?del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Unknown offset should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_unknown_offset_single_position() {
        // Even a single unknown offset should cause early return
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10-?del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Single unknown offset should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_utr_before_tx_start_returns_unchanged() {
        // c.-215 with a small UTR should not error - return unchanged
        // NM_001234.1 has cds_start=5, so c.-215 maps to 5 + (-215) - 1 = -211
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.-215_-214del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "UTR before transcript start should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_no_cds_returns_unchanged() {
        // An NR_ transcript with c. coordinates should not error
        let mut provider = MockProvider::new();
        use crate::reference::transcript::{Exon, Transcript};
        provider.add_transcript(Transcript {
            id: "NR_001566.1".to_string(),
            gene_symbol: Some("NCRNA".to_string()),
            strand: crate::reference::transcript::Strand::Plus,
            sequence: Some("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 51)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: Default::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });
        let normalizer = Normalizer::new(provider);

        // c. variant on a non-coding transcript (no CDS)
        let variant = parse_hgvs("NR_001566.1:c.10del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "No CDS should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_tx_intronic() {
        // n. intronic variants should normalize via genomic space
        // Build a non-coding transcript with genomic coords and intronic positions
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};

        let mut provider = MockProvider::new();

        // Create transcript: 2 exons with an intron in between
        // Exon 1: tx 1-100, genomic 1000-1099
        // Intron: genomic 1100-1199
        // Exon 2: tx 101-200, genomic 1200-1299
        // Sequence in the intron around position 1100+: AAAA... (for shifting test)
        let tx_sequence = "A".repeat(200);

        provider.add_transcript(Transcript {
            id: "NR_038982.1".to_string(),
            gene_symbol: Some("NCRNA_TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_sequence),
            cds_start: None,
            cds_end: None,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 1200, 1299),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1299),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });

        // Add genomic sequence for chr1 around positions 1000-1299
        // Make the intron region (1100-1199) be "AGCT" repeated to test shifting
        let mut genomic = String::new();
        for _ in 0..325 {
            genomic.push_str("AGCT");
        }
        provider.add_genomic_sequence("chr1", genomic);

        let normalizer = Normalizer::new(provider);

        // n.100+4del - intronic deletion in a non-coding transcript
        let variant = parse_hgvs("NR_038982.1:n.100+4del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "n. intronic normalization should succeed, got: {:?}",
            result.err()
        );
        let output = format!("{}", result.unwrap());
        assert!(
            output.contains('+') || output.contains('-'),
            "Normalized intronic n. variant should retain intronic notation, got: {}",
            output
        );
    }
}
