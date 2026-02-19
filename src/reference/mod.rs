//! Reference data abstraction
//!
//! Provides traits and implementations for accessing reference sequence data.

pub mod fasta;
pub mod loader;
pub mod mock;
pub mod multi_fasta;
pub mod protein;
pub mod provider;
pub mod transcript;

pub use fasta::FastaProvider;
pub use loader::{detect_genome_build, load_gff3, load_gtf, TranscriptDb};
pub use mock::MockProvider;
pub use multi_fasta::MultiFastaProvider;
pub use protein::ProteinCache;
pub use provider::ReferenceProvider;
pub use transcript::{Exon, Strand, Transcript};
