//! HGVS types and parser
//!
//! This module contains the complete type system for representing HGVS variants
//! and a nom-based parser for parsing HGVS strings.

pub mod edit;
pub mod interval;
pub mod location;
pub mod parser;
pub mod uncertainty;
pub mod validation;
pub mod variant;
pub mod version;

// Re-export commonly used types
pub use edit::{FrameshiftTer, NaEdit, ProteinEdit, RepeatUnit};
pub use interval::{CdsInterval, GenomeInterval, ProtInterval, TxInterval};
pub use location::{
    AaCode, CdsPos, GenomePos, ProtPos, ProteinRenderStyle, SpecialPosition, TerStyle, TxPos,
};
pub use uncertainty::Mu;
pub use validation::{ParseConfig, ValidationLevel};
pub use variant::{AllelePhase, AlleleVariant, CoordinateAxis, HgvsVariant};
pub use version::HgvsVersion;
