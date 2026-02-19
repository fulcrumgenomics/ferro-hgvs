//! Data loading and coordinate mapping.
//!
//! This module provides functionality for loading transcript alignment data
//! (cdot format) and mapping coordinates between different coordinate systems.
//!
//! # Coordinate Mapping
//!
//! The main entry point is [`CoordinateMapper`] which provides methods for
//! converting between CDS (c.) and genomic (g.) coordinates.
//!
//! # Example
//!
//! ```no_run
//! use ferro_hgvs::data::{CdotMapper, CoordinateMapper};
//! use ferro_hgvs::hgvs::location::CdsPos;
//!
//! // Load transcript data
//! let cdot = CdotMapper::from_json_file("transcripts.json").unwrap();
//! let mapper = CoordinateMapper::new(cdot);
//!
//! // Convert c.100 to genomic coordinates
//! let cds_pos = CdsPos::new(100);
//! let result = mapper.cds_to_genome("NM_000088.3", &cds_pos).unwrap();
//! println!("Genomic position: {}", result.variant.base);
//! ```

pub mod cdot;
pub mod mapping;
pub mod projection;

pub use cdot::{CdotFile, CdotMapper, CdotTranscript, CdsPosition, Exon};
pub use mapping::{CoordinateMapper, MappingInfo, MappingResult};
pub use projection::{ManeStatus, ProjectionResult, Projector, TranscriptProjection};
