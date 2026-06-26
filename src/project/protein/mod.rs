//! Protein consequence prediction for CDS variants.
//!
//! This module predicts the amino acid level consequence of a variant using the
//! HGVS protein nomenclature recommendations:
//! <https://hgvs-nomenclature.org/stable/recommendations/protein/>
//!
//! # Submodules
//!
//! * [`substitution`] — substitution (missense / synonymous / nonsense) prediction.
//! * [`indel`]         — deletion / insertion / duplication / delins / inversion prediction.
//! * [`helpers`]       — shared pure helpers (CDS building, translation, diff-finding).

mod helpers;
mod indel;
mod substitution;

// Re-export the public API used by `projector.rs`.
pub(crate) use helpers::{
    build_initiator_unknown, build_whole_protein_unknown, cds_has_recognized_start,
    cds_has_valid_start, edit_reaches_initiation_codon, read_cds_start_codon,
    whole_exon_deletion_span, RefProteinBundle,
};
pub(crate) use indel::predict_indel_protein;
// `read_ref_codon` / `translate` are the position-level CDS primitives reused by
// the service effect handler to resolve real amino-acid residues (issue #806);
// `read_ref_codon` is also the seam #498 (full c.→p.) inherits.
pub(crate) use substitution::{predict_substitution_protein, read_ref_codon, translate};
