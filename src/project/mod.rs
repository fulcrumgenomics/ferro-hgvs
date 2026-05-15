//! Variant-level projection: g. variants → c./p. equivalents on a target transcript.

mod accession;
mod edit;
mod projector;
mod protein;
mod result;

pub use projector::VariantProjector;
pub use result::VariantProjection;
