//! Variant-level projection: g. variants → c./p. equivalents on a target transcript.

mod accession;
pub mod edit;
mod projector;
mod protein;
mod result;
pub(crate) mod rna;
mod transcript_axis;

pub use projector::VariantProjector;
pub use result::VariantProjection;
