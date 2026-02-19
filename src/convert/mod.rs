//! Coordinate conversion
//!
//! Provides conversion between different HGVS coordinate systems:
//! - Genomic (g.) ↔ Transcript
//! - CDS (c.) ↔ Transcript (n.)
//! - CDS (c.) → Protein (p.)

pub mod coding;
pub mod genomic;
pub mod mapper;
pub mod noncoding;
pub mod protein;

pub use mapper::CoordinateMapper;
