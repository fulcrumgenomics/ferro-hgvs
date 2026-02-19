//! MAVE-HGVS support for MaveDB variant notation.
//!
//! This module provides context-aware parsing for MAVE-HGVS short-form notation
//! used by MaveDB (Multiplex Assay of Variant Effect database).
//!
//! # MAVE-HGVS Format
//!
//! MAVE-HGVS is a simplified form of HGVS that often omits reference accessions:
//! - `p.Leu11Pro` instead of `NP_000509:p.Leu11Pro`
//! - `c.32A>G` instead of `NM_000518:c.32A>G`
//!
//! This short form is valid within the context of a specific score set that
//! defines the target sequence.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::mave::{MaveContext, parse_mave_hgvs};
//!
//! // Create a context with the target sequence
//! let context = MaveContext::new()
//!     .with_protein_accession("NP_000509.1")
//!     .with_coding_accession("NM_000518.5")
//!     .with_gene_symbol("HBB");
//!
//! // Parse a protein variant
//! let variant = parse_mave_hgvs("p.Glu6Val", &context).unwrap();
//! assert_eq!(variant.to_string(), "NP_000509.1:p.Glu6Val");
//!
//! // Parse a coding variant
//! let variant = parse_mave_hgvs("c.20A>T", &context).unwrap();
//! assert_eq!(variant.to_string(), "NM_000518.5:c.20A>T");
//! ```
//!
//! # References
//!
//! - [MaveDB](https://www.mavedb.org/)
//! - [MAVE-HGVS Specification](https://mavedb.org/docs/mavehgvs/)

mod context;
mod parser;

pub use context::MaveContext;
pub use parser::{is_mave_short_form, parse_mave_hgvs, parse_mave_hgvs_lenient, MaveParseError};
