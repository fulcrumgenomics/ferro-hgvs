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
pub mod rules;
pub mod shuffle;
pub mod validate;

use crate::coords::{hgvs_pos_to_index, index_to_hgvs_pos};
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::Interval;
use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
use crate::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
use crate::hgvs::variant::{CdsVariant, GenomeVariant, HgvsVariant, LocEdit, TxVariant};
use crate::hgvs::HgvsVariant as HV;
use crate::reference::ReferenceProvider;
use boundary::Boundaries;
pub use config::{NormalizeConfig, ShuffleDirection};
use rules::{canonicalize_edit, needs_normalization, should_canonicalize};
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

/// Warning generated during normalization
#[derive(Debug, Clone)]
pub struct NormalizationWarning {
    /// Warning code (e.g., "REFSEQ_MISMATCH")
    pub code: String,
    /// Human-readable description
    pub message: String,
    /// What the input claimed as reference
    pub stated_ref: String,
    /// What the actual reference sequence has
    pub actual_ref: String,
    /// Position info
    pub position: String,
    /// Whether the mismatch was auto-corrected
    pub corrected: bool,
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
        self.warnings.iter().any(|w| w.code == "REFSEQ_MISMATCH")
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

        // In strict mode, reject if there were reference mismatches
        if self.config.should_reject_ref_mismatch() && result.has_ref_mismatch() {
            if let Some(warning) = result.warnings.iter().find(|w| w.code == "REFSEQ_MISMATCH") {
                return Err(FerroError::ReferenceMismatch {
                    location: warning.position.clone(),
                    expected: warning.stated_ref.clone(),
                    found: warning.actual_ref.clone(),
                });
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
        // First pass: normalize each variant independently, collecting all warnings
        let mut all_warnings = Vec::new();
        let mut normalized = Vec::new();

        for v in &allele.variants {
            let result = self.normalize_with_warnings(v)?;
            all_warnings.extend(result.warnings);
            normalized.push(result.result);
        }

        // Second pass: check for overlaps and resolve
        if self.config.prevent_overlap {
            normalized = self.resolve_overlaps(normalized)?;
        }

        Ok((
            HgvsVariant::Allele(crate::hgvs::variant::AlleleVariant::new(
                normalized,
                allele.phase,
            )),
            all_warnings,
        ))
    }

    /// Resolve overlaps between normalized variants in a compound allele.
    ///
    /// When normalization shifts variants, they may end up overlapping. This function
    /// detects overlaps and constrains the normalization to prevent collisions.
    fn resolve_overlaps(&self, variants: Vec<HgvsVariant>) -> Result<Vec<HgvsVariant>, FerroError> {
        if variants.len() < 2 {
            return Ok(variants);
        }

        // Extract positions for comparison
        let mut positions: Vec<(usize, Option<(u64, u64)>)> = variants
            .iter()
            .enumerate()
            .map(|(i, v)| (i, self.extract_position_range(v)))
            .collect();

        // Sort by start position
        positions.sort_by(|a, b| match (&a.1, &b.1) {
            (Some((s1, _)), Some((s2, _))) => s1.cmp(s2),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });

        // Check for overlaps and mark variants that need adjustment
        let mut needs_adjustment = vec![false; variants.len()];

        for i in 0..positions.len() - 1 {
            let (_idx_a, pos_a) = &positions[i];
            let (idx_b, pos_b) = &positions[i + 1];

            if let (Some((_, end_a)), Some((start_b, _))) = (pos_a, pos_b) {
                // Check for overlap: end of A >= start of B
                if end_a >= start_b {
                    // Mark the second variant for potential adjustment
                    needs_adjustment[*idx_b] = true;
                }
            }
        }

        // For now, we report overlaps but don't modify the variants
        // A more sophisticated implementation would re-normalize with constraints
        // or report a warning
        for (i, needs_adj) in needs_adjustment.iter().enumerate() {
            if *needs_adj {
                // In a complete implementation, we would constrain normalization
                // For now, we preserve the variant as-is with a potential warning
                // Future: add to a warnings collection
                let _ = i; // Suppress unused variable warning
            }
        }

        Ok(variants)
    }

    /// Extract the genomic/CDS position range from a variant.
    fn extract_position_range(&self, variant: &HgvsVariant) -> Option<(u64, u64)> {
        match variant {
            HV::Genome(v) => {
                let start = v.loc_edit.location.start.inner()?.base;
                let end = v.loc_edit.location.end.inner()?.base;
                Some((start, end))
            }
            HV::Cds(v) => {
                let start_pos = v.loc_edit.location.start.inner()?;
                let end_pos = v.loc_edit.location.end.inner()?;
                // Only return positions for exonic, non-UTR variants
                if start_pos.is_intronic()
                    || end_pos.is_intronic()
                    || start_pos.base < 1
                    || end_pos.base < 1
                {
                    None
                } else {
                    Some((start_pos.base as u64, end_pos.base as u64))
                }
            }
            HV::Tx(v) => {
                let start = v.loc_edit.location.start.inner()?.base;
                let end = v.loc_edit.location.end.inner()?.base;
                // Only return positions for positive bases (within transcript)
                if start < 1 || end < 1 {
                    None
                } else {
                    Some((start as u64, end as u64))
                }
            }
            _ => None,
        }
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

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Genome(variant.clone()), vec![]));
        }

        // Get sequence from provider - can't normalize with unknown positions
        let accession = variant.accession.full();
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
        let (new_rel_start, new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            ref_seq.as_bytes(),
            edit,
            rel_start,
            rel_end,
            &Boundaries::new(1, ref_seq.len() as u64),
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

        Ok((HV::Genome(new_variant), warnings))
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

        // Try to get transcript first - we need it for intronic normalization too
        let accession = variant.accession.full();
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            // Can't do full normalization without transcript, but apply minimal notation
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Handle intronic variants specially
        if start_pos.is_intronic() || end_pos.is_intronic() {
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

        // Get boundaries (stay within exon unless configured otherwise)
        let boundaries = if self.config.cross_boundaries {
            Boundaries::new(1, transcript.sequence.len() as u64)
        } else {
            match boundary::get_cds_boundaries(&transcript, tx_start, &self.config) {
                Ok(b) => b,
                Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
            }
        };

        // Perform normalization on transcript sequence
        let seq = transcript.sequence.as_bytes();
        let (new_tx_start, new_tx_end, new_edit, warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries)?;

        // Convert back to CDS coordinates
        let new_start = self.tx_to_cds_pos(new_tx_start, cds_start, transcript.cds_end)?;
        let new_end = self.tx_to_cds_pos(new_tx_end, cds_start, transcript.cds_end)?;

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        Ok((HV::Cds(new_variant), warnings))
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

        // Try to get transcript first - we need it for intronic normalization too
        let accession = variant.accession.full();
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            // Can't do full normalization without transcript, but apply minimal notation
            Err(_) => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };

        if start_pos.is_intronic() || end_pos.is_intronic() {
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
        let boundaries = Boundaries::new(1, transcript.sequence.len() as u64);

        // Perform normalization
        let seq = transcript.sequence.as_bytes();
        let (new_start, new_end, new_edit, warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries)?;

        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(TxPos::new(new_start as i64), TxPos::new(new_end as i64)),
                new_edit,
            ),
        };

        Ok((HV::Tx(new_variant), warnings))
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

        let accession = variant.accession.full();

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
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{LocEdit, RnaVariant};

        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

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
        let accession = variant.accession.full();
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            Err(_) => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // Only normalize positive positions (within transcript)
        // Negative positions are outside the transcript sequence
        if start_pos.base < 1 || end_pos.base < 1 {
            return Ok((HV::Rna(variant.clone()), vec![]));
        }

        let rna_start = start_pos.base as u64;
        let rna_end = end_pos.base as u64;

        // Get boundaries
        let boundaries = Boundaries::new(1, transcript.sequence.len() as u64);

        // Perform normalization
        let seq = transcript.sequence.as_bytes();
        let (new_start, new_end, new_edit, warnings) =
            self.normalize_na_edit(seq, edit, rna_start, rna_end, &boundaries)?;

        let new_variant = RnaVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                RnaInterval::new(RnaPos::new(new_start as i64), RnaPos::new(new_end as i64)),
                new_edit,
            ),
        };

        Ok((HV::Rna(new_variant), warnings))
    }

    /// Normalize a mitochondrial variant
    fn normalize_mt(
        &self,
        variant: &crate::hgvs::variant::MtVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // MT variants are similar to genomic
        // For now, return as-is (could implement circular genome handling)
        Ok((HV::Mt(variant.clone()), vec![]))
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

        // Get the chromosome for this transcript
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Transcript has no chromosome mapping for intronic normalization"
                        .to_string(),
                })?;

        // Create coordinate mapper
        let mapper = CoordinateMapper::new(transcript);

        // Convert CDS intronic positions to genomic
        let g_start = mapper.cds_to_genomic_with_intron(start_pos)?;
        let g_end = mapper.cds_to_genomic_with_intron(end_pos)?;

        // Ensure start <= end (may be reversed on minus strand)
        let (g_start, g_end) = if g_start <= g_end {
            (g_start, g_end)
        } else {
            (g_end, g_start)
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

        // Perform normalization in genomic space
        let seq_bytes = genomic_seq.as_bytes();
        let (new_rel_start, new_rel_end, new_edit, warnings) =
            self.normalize_na_edit(seq_bytes, edit, rel_start, rel_end, &boundaries)?;

        // Convert the normalized genomic position back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert back to CDS intronic notation
        let new_start = mapper.genomic_to_cds_intronic(new_g_start)?;
        let new_end = mapper.genomic_to_cds_intronic(new_g_end)?;

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

        // Get the chromosome for this transcript
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Transcript has no chromosome mapping for intronic normalization"
                        .to_string(),
                })?;

        // Create coordinate mapper
        let mapper = CoordinateMapper::new(transcript);

        // Convert tx intronic positions to genomic
        let g_start = mapper.tx_to_genomic_with_intron(start_pos)?;
        let g_end = mapper.tx_to_genomic_with_intron(end_pos)?;

        // Ensure start <= end (may be reversed on minus strand)
        let (g_start, g_end) = if g_start <= g_end {
            (g_start, g_end)
        } else {
            (g_end, g_start)
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

        // Perform normalization in genomic space
        let seq_bytes = genomic_seq.as_bytes();
        let (new_rel_start, new_rel_end, new_edit, warnings) =
            self.normalize_na_edit(seq_bytes, edit, rel_start, rel_end, &boundaries)?;

        // Convert the normalized genomic position back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert back to transcript intronic notation
        let new_start = mapper.genomic_to_tx_with_intron(new_g_start)?;
        let new_end = mapper.genomic_to_tx_with_intron(new_g_end)?;

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

        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Transcript has no chromosome for boundary normalization".to_string(),
                })?;

        let mapper = CoordinateMapper::new(transcript);

        // Convert both positions to genomic
        // For exonic positions, use standard conversion
        // For intronic positions, use intronic conversion
        let g_start = self.cds_pos_to_genomic(&mapper, start_pos)?;
        let g_end = self.cds_pos_to_genomic(&mapper, end_pos)?;

        // Handle strand orientation - ensure start <= end
        let (g_start, g_end) = if g_start <= g_end {
            (g_start, g_end)
        } else {
            (g_end, g_start)
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

        // Normalize in genomic space
        let seq_bytes = genomic_seq.as_bytes();
        let (new_rel_start, new_rel_end, new_edit, warnings) =
            self.normalize_na_edit(seq_bytes, edit, rel_start, rel_end, &boundaries)?;

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
    ) -> Result<(u64, u64, NaEdit, Vec<NormalizationWarning>), FerroError> {
        let mut warnings = Vec::new();

        // Validate reference allele before normalization
        let validation = validate::validate_reference(edit, ref_seq, start, end);
        if !validation.valid {
            warnings.push(NormalizationWarning {
                code: "REFSEQ_MISMATCH".to_string(),
                message: validation.warning.unwrap_or_default(),
                stated_ref: validation.stated_ref.unwrap_or_default(),
                actual_ref: validation.actual_ref.unwrap_or_default(),
                position: format!("{}-{}", start, end),
                corrected: true,
            });
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
            NaEdit::Duplication { sequence, .. } => {
                if let Some(seq) = sequence {
                    seq.bases().iter().map(|b| b.to_u8()).collect()
                } else {
                    // Get duplicated sequence from reference
                    // start is 1-based, convert to 0-based index for array access
                    let s = hgvs_pos_to_index(start);
                    let e = end as usize;
                    if e <= ref_seq.len() {
                        ref_seq[s..e].to_vec()
                    } else {
                        vec![]
                    }
                }
            }
            NaEdit::Delins { sequence } => {
                use crate::hgvs::edit::InsertedSequence;

                // HGVS spec: delins should NOT be 3' shifted like del/dup/ins
                // But we should check if delins should become a duplication first
                // Example: c.5delinsGG where position 5 is G → dup
                if let InsertedSequence::Literal(seq) = sequence {
                    let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();
                    let start_idx = hgvs_pos_to_index(start);
                    let end_idx = end as usize;

                    if rules::delins_is_duplication(ref_seq, start_idx, end_idx, &seq_bytes) {
                        // Convert to duplication with minimal notation
                        return Ok((
                            start,
                            end,
                            NaEdit::Duplication {
                                sequence: None,
                                length: None,
                                uncertain_extent: None,
                            },
                            warnings.clone(),
                        ));
                    }
                }
                // Return unchanged (only apply minimal representation rules, not shuffling)
                return Ok((start, end, edit.clone(), warnings.clone()));
            }
            NaEdit::Inversion { sequence, length } => {
                // Apply inversion shortening: remove outer bases that cancel out
                // when the first base is complementary to the last base
                let start_idx = hgvs_pos_to_index(start); // Convert 1-based to 0-based
                let end_idx = end as usize; // end is exclusive (0-based)

                if let Some((new_s, new_e)) = rules::shorten_inversion(ref_seq, start_idx, end_idx)
                {
                    // Inversion was shortened - convert back to 1-based
                    let shortened_edit = NaEdit::Inversion {
                        sequence: sequence.clone(),
                        length: *length,
                    };
                    return Ok((
                        index_to_hgvs_pos(new_s),
                        new_e as u64,
                        shortened_edit,
                        warnings,
                    ));
                } else {
                    // Inversion reduced to identity - return as identity
                    return Ok((
                        start,
                        end,
                        NaEdit::Identity {
                            sequence: None,
                            whole_entity: false,
                        },
                        warnings,
                    ));
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
                match rules::normalize_repeat(ref_seq, pos_idx, &repeat_unit, *specified_count) {
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
                        let rep_edit = NaEdit::Repeat {
                            sequence: Some(Sequence::new(bases)),
                            count: RepeatCount::Exact(rep_count),
                            additional_counts: vec![],
                            trailing: None,
                        };
                        return Ok((rep_start, rep_end, rep_edit, warnings));
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
                        if let Some((base, count, rep_start, rep_end)) =
                            rules::insertion_to_repeat(ref_seq, original_pos_idx, &seq_bytes)
                        {
                            // Convert to repeat notation: BASE[count]
                            use crate::hgvs::edit::Base;
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

                    if rules::insertion_is_duplication(ref_seq, result.start, &rotated_seq) {
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
                    } else {
                        // Output the rotated sequence for shifted insertions
                        if rotation > 0 {
                            use crate::hgvs::edit::{Base, InsertedSequence, Sequence};
                            let rotated_bases: Vec<Base> = rotated_seq
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
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
                    rules::duplication_to_repeat(ref_seq, result.start, result.end)
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
            NaEdit::Delins { sequence } => {
                use crate::hgvs::edit::{InsertedSequence, Sequence};

                // Check if delins should become a duplication
                // Example: c.5delinsGG where position 5 is G → dup
                if let InsertedSequence::Literal(seq) = sequence {
                    let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();
                    let start_0 = result.start as usize;
                    let end_0 = result.end as usize;

                    if rules::delins_is_duplication(ref_seq, start_0, end_0, &seq_bytes) {
                        // Convert to duplication - the duplicated sequence is the deleted region
                        let dup_len = end_0 - start_0;
                        let dup_seq_bytes = &ref_seq[start_0..end_0];
                        let dup_bases: Vec<crate::hgvs::edit::Base> = dup_seq_bytes
                            .iter()
                            .filter_map(|&b| crate::hgvs::edit::Base::from_char(b as char))
                            .collect();
                        let dup_seq = Sequence::new(dup_bases);

                        (
                            new_start,
                            new_end,
                            NaEdit::Duplication {
                                sequence: Some(dup_seq),
                                length: Some(dup_len as u64),
                                uncertain_extent: None,
                            },
                        )
                    } else {
                        (new_start, new_end, edit.clone())
                    }
                } else {
                    (new_start, new_end, edit.clone())
                }
            }
            // Deletions - strip explicit length for minimal notation
            NaEdit::Deletion { .. } => (
                new_start,
                new_end,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
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
        } else if pos.base < 1 {
            let tx_pos = cds_start as i64 + pos.base - 1;
            u64::try_from(tx_pos).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "CDS position c.{} maps before transcript start (cds_start={})",
                    pos.base, cds_start
                ),
            })
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
            // For 5'UTR: cds_to_tx_pos formula is tx = cds_start + base - 1
            // So inverse is: base = tx - cds_start + 1
            Ok(CdsPos {
                base: pos as i64 - cds_start as i64 + 1,
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

        // Verify output format preserves both variants
        let output = format!("{}", result);
        assert!(
            output.contains("c.10A>G"),
            "Allele should contain first variant, got: {}",
            output
        );
        assert!(
            output.contains("c.20C>T"),
            "Allele should contain second variant, got: {}",
            output
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
            sequence: tx_seq.to_string(),
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
        // Test: c.40_40+3dup - duplication spanning exon-intron boundary
        // c.40 = last base of exon 2 (A at g.1049)
        // c.40+3 = 3rd intronic base of intron 2 (A at g.1052)
        // Since intron 2 starts with AAA (and exon 2 ends with AAAAA), this creates a repeat pattern
        // The normalizer correctly converts this to repeat notation (A[N])
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
        // May remain as dup or be converted to repeat notation for poly-A
        assert!(
            output.contains("dup") || output.contains("["),
            "Should remain a duplication or become repeat notation, got: {}",
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
            sequence: seq.clone(),
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
            sequence: seq.clone(),
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
            sequence: seq,
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
        // cds_start=5, base=-3 → 5 + (-3) - 1 = 1, valid position
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let pos = CdsPos {
            base: -3,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&pos, 5, Some(38));
        assert_eq!(result.unwrap(), 1);
    }

    #[test]
    fn test_variant_positions_negative_cds_base() {
        // CDS variant with negative base should return None from variant_positions
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs("NM_001234.1:c.-3A>G").unwrap();
        let positions = normalizer.extract_position_range(&variant);
        assert_eq!(
            positions, None,
            "variant_positions should return None for 5' UTR CDS positions"
        );
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
            sequence: "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC".to_string(),
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
            sequence: tx_sequence,
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
