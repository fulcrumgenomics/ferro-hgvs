//! Reference data abstraction
//!
//! Provides traits and implementations for accessing reference sequence data.

pub mod annotation;
pub mod authoritative;
pub mod derived_placement;
pub mod derived_tx_structure;
pub mod fasta;
pub mod legacy_selector;
pub mod loader;
pub mod mock;
pub mod multi_fasta;
pub mod ng_hosted_transcripts;
pub mod ng_placement_builder;
pub mod protein;
pub mod provider;
pub mod transcript;
pub mod validate;

pub use annotation::{load_annotations, AnnotationFormat, LoaderConfig, LoaderReport};
pub use fasta::FastaProvider;
pub use loader::{detect_genome_build, TranscriptDb};
pub use mock::{JsonProvider, MockProvider};
pub use multi_fasta::MultiFastaProvider;
pub use protein::ProteinCache;
pub use provider::{GenomicPlacement, ReferenceProvider};
pub use transcript::{Exon, Strand, Transcript};
