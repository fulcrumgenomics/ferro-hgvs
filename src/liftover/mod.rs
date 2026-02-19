//! Cross-genome liftover for coordinate conversion between genome builds.
//!
//! This module provides functionality to convert genomic coordinates between
//! GRCh37 (hg19) and GRCh38 (hg38) using UCSC chain files.
//!
//! # Overview
//!
//! Liftover is the process of converting coordinates from one genome assembly
//! to another. This is necessary because different genome builds have slightly
//! different coordinate systems due to assembly improvements.
//!
//! # Example
//!
//! ```no_run
//! use ferro_hgvs::liftover::{ChainFile, Liftover};
//! use ferro_hgvs::reference::transcript::GenomeBuild;
//!
//! // Load chain files
//! let forward = ChainFile::from_file("hg19ToHg38.over.chain.gz").unwrap();
//! let reverse = ChainFile::from_file("hg38ToHg19.over.chain.gz").unwrap();
//! let liftover = Liftover::new(forward, reverse);
//!
//! // Lift a position from GRCh37 to GRCh38
//! let result = liftover.lift(
//!     GenomeBuild::GRCh37,
//!     GenomeBuild::GRCh38,
//!     "chr1",
//!     12345,
//! ).unwrap();
//!
//! println!("Lifted to {}:{}", result.target_contig, result.target_pos);
//! ```
//!
//! # Chain File Format
//!
//! UCSC chain files define alignment blocks between genome assemblies.
//! Chain files can be downloaded from:
//! - <https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz>
//! - <https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz>

pub mod aliases;
pub mod chain;
pub mod lift;

pub use aliases::ContigAliases;
pub use chain::{Chain, ChainBlock, ChainFile};
pub use lift::{IntervalLiftoverResult, Liftover, LiftoverResult, LiftoverStatus};
