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

use crate::FerroError;

/// Reject an on-disk artifact whose schema version is newer than this build
/// supports. Absent versions deserialize to 0 and are accepted as legacy.
///
/// Shared by the derived-artifact loaders (`derived_placement`,
/// `derived_tx_structure`) — it lives at the `reference` module root rather
/// than in either sibling because neither owns the other, and a single
/// definition is what keeps the operator-facing wording from drifting apart
/// between the two artifacts (#1001).
///
/// * `artifact` — the manifest key naming the artifact (e.g.
///   `"derived_transcript_placements"`), so the message points at what to fix.
/// * `found` / `max` — the version read from disk and the maximum this build reads.
/// * `path` — the artifact's path, echoed into the message.
pub fn check_artifact_schema_version(
    artifact: &str,
    found: u32,
    max: u32,
    path: &std::path::Path,
) -> Result<(), FerroError> {
    if found > max {
        return Err(FerroError::Io {
            msg: format!(
                "{} artifact at {} has schema_version {}, \
                 which is newer than this build supports (maximum {}). Upgrade ferro, or \
                 re-derive the artifact.",
                artifact,
                path.display(),
                found,
                max
            ),
        });
    }
    Ok(())
}

pub use annotation::{load_annotations, AnnotationFormat, LoaderConfig, LoaderReport};
pub use fasta::FastaProvider;
pub use loader::{detect_genome_build, TranscriptDb};
pub use mock::{JsonProvider, MockProvider};
pub use multi_fasta::MultiFastaProvider;
pub use protein::ProteinCache;
pub use provider::{GenomicPlacement, ReferenceProvider};
pub use transcript::{Exon, Strand, Transcript};
