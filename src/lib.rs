// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! ferro-hgvs: HGVS variant normalizer
//!
//! Part of the ferro bioinformatics toolkit.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::{parse_hgvs, Normalizer, JsonProvider};
//!
//! // Parse an HGVS variant string
//! let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
//!
//! // Create a normalizer with test data
//! let provider = JsonProvider::with_test_data();
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
pub mod conformance;
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
pub mod parallel;
/// Performance-comparison table types + rendering (used by the
/// `generate_perf_tables` example and its tests).
#[cfg(feature = "dev")]
pub mod perf_table;
pub mod prepare;
pub mod project;
#[cfg(feature = "python")]
pub mod python;
pub mod python_helpers;
pub mod reference;
pub mod rsid;
pub mod sequence;
#[cfg(feature = "web-service")]
pub mod service;
pub mod spdi;
/// Tool-support matrix types + table generation (used by the
/// `generate_tool_support_tables` example and its tests).
#[cfg(feature = "dev")]
pub mod tool_support;
pub mod vcf;

// Re-export commonly used types
pub use error::FerroError;
pub use hgvs::location::{AaCode, ProteinRenderStyle, TerStyle};
pub use hgvs::parser::{parse_hgvs, parse_hgvs_fast};
pub use hgvs::variant::{CoordinateAxis, HgvsVariant};
pub use normalize::from_sequences::{
    from_sequences, from_sequences_detailed, DerivedDescription, FromSequencesOptions,
};
pub use normalize::sequence_pair::SequencePair;
pub use normalize::{NormalizeConfig, Normalizer};
// `ShuffleDirection` is re-exported for ferro's own integration tests only —
// `tests/` is an external crate, so the 3'/5' differential oracle cannot reach a
// `pub(crate)` type. It is `#[doc(hidden)]` at its definition and no ferro entry
// point accepts a direction from a caller; see the type's docs and `README.md`
// rule 6.
#[doc(hidden)]
pub use normalize::ShuffleDirection;
pub use project::{VariantProjection, VariantProjector};
pub use reference::{JsonProvider, MockProvider, MultiFastaProvider, ReferenceProvider};
pub use spdi::{
    hgvs_to_spdi_simple, parse_spdi, spdi_to_hgvs, spdi_to_hgvs_with_ref, ConversionError,
    SpdiVariant,
};

// Re-export coordinate types for type-safe position handling
pub use coords::{
    cdot_genomic_to_closed, cdot_tx_coords, hgvs_pos_to_index, hgvs_pos_to_index_checked,
    hgvs_to_spdi_pos, hgvs_to_spdi_pos_checked, index_to_hgvs_pos, spdi_to_hgvs_pos,
    OneBasedInterval, OneBasedPos, ZeroBasedInterval, ZeroBasedPos,
};

/// Result type alias for ferro-hgvs operations
pub type Result<T> = std::result::Result<T, FerroError>;
