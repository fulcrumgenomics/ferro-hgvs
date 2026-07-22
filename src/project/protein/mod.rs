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
//! * [`identity`]      — describing an all-silent change by the codons it rewrote.
//! * [`helpers`]       — shared pure helpers (CDS building, translation, diff-finding).

mod helpers;
mod identity;
mod indel;
mod substitution;

// Re-export the public API used by `projector.rs`.
pub(crate) use helpers::{
    build_initiator_unknown, build_whole_protein_unknown, cds_has_recognized_start,
    cds_has_valid_start, edit_reaches_initiation_codon, edit_spans_cds_into_3utr,
    read_cds_start_codon, whole_exon_deletion_span, RefProteinBundle,
};
// The one place an all-silent change is turned into a description (#1099), so
// the codon-level and residue-level cis combiners in `projector.rs` and the
// indel path cannot drift apart in how they name the codons they rewrote.
pub(crate) use identity::{codon_residues, group_consecutive_by, render_silent_identity};
pub(crate) use indel::{
    combined_cis_residue_changes, predict_indel_protein, predict_stop_region_extension,
    try_project_cis_combined_inframe, CisCombined,
};
pub(crate) use substitution::predict_substitution_protein;
// `apply_substitution` / `read_ref_codon` / `translate` are the position-level
// CDS primitives reused by the projector's codon-level in-cis substitution
// combiner (#1076) and by the service effect handler (#806); `read_ref_codon`
// is also the seam #498 (full c.→p.) inherits.
pub(crate) use substitution::{apply_substitution, read_ref_codon, translate};
