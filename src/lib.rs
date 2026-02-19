// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! ferro-hgvs: HGVS variant normalizer
//!
//! Part of the ferro bioinformatics toolkit.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::{parse_hgvs, Normalizer, MockProvider};
//!
//! // Parse an HGVS variant string
//! let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
//!
//! // Create a normalizer with test data
//! let provider = MockProvider::with_test_data();
//! let normalizer = Normalizer::new(provider);
//!
//! // Normalize the variant
//! let normalized = normalizer.normalize(&variant).unwrap();
//! println!("Normalized: {}", normalized);
//! ```

pub mod backtranslate;
pub mod batch;
#[cfg(feature = "benchmark")]
pub mod benchmark;
pub mod cache;
pub mod check;
pub mod cli;
pub mod clinvar;
pub mod commands;
pub mod config;
pub mod convert;
pub mod coords;
pub mod data;
pub mod diagnostic;
pub mod effect;
pub mod equivalence;
pub mod error;
pub mod error_handling;
pub mod extractor;
pub mod hgvs;
pub mod legacy;
pub mod liftover;
pub mod mave;
pub mod normalize;
#[cfg(feature = "parallel")]
pub mod parallel;
pub mod prepare;
#[cfg(feature = "python")]
pub mod python;
pub mod python_helpers;
pub mod reference;
pub mod rsid;
#[cfg(feature = "web-service")]
pub mod service;
pub mod spdi;
pub mod vcf;

// Re-export commonly used types
pub use error::FerroError;
pub use hgvs::parser::{parse_hgvs, parse_hgvs_fast};
pub use hgvs::variant::HgvsVariant;
pub use normalize::{NormalizeConfig, Normalizer, ShuffleDirection};
pub use reference::{MockProvider, MultiFastaProvider, ReferenceProvider};
pub use spdi::{hgvs_to_spdi_simple, parse_spdi, spdi_to_hgvs, ConversionError, SpdiVariant};

// Re-export coordinate types for type-safe position handling
pub use coords::{
    cdot_genomic_to_closed, cdot_tx_coords, hgvs_pos_to_index, hgvs_to_spdi_pos, index_to_hgvs_pos,
    spdi_to_hgvs_pos, OneBasedInterval, OneBasedPos, ZeroBasedInterval, ZeroBasedPos,
};

/// Result type alias for ferro-hgvs operations
pub type Result<T> = std::result::Result<T, FerroError>;
