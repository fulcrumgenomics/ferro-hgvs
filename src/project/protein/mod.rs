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
pub(crate) use helpers::RefProteinBundle;
pub(crate) use indel::predict_indel_protein;
pub(crate) use substitution::predict_substitution_protein;
