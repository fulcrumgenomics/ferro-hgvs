//! Variant-level projection: g. variants → c./p. equivalents on a target transcript.

pub(crate) mod accession;
pub(crate) mod codon_exception;
pub mod edit;
mod mutalyzer;
mod projector;
pub(crate) mod protein;
mod result;
pub(crate) mod rna;
mod transcript_axis;

pub use mutalyzer::{MutalyzerResult, MutalyzerWarning};
pub use projector::VariantProjector;
pub use result::{AxisDeclineReasons, VariantProjection};
