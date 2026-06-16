//! Tests for issue #314 — `MmapFastaProvider::get_transcript` carries the
//! same permissive `starts_with(base_id)` fallback that `MockProvider`
//! had before #312. A lookup for a chromosome accession (e.g. `chr1`)
//! resolves to any transcript whose id literally starts with that
//! string (e.g. `chr1-gene.1`).
//!
//! These tests pin the corrected lookup contract:
//!
//! - the unversioned-prefix fallback must match a stored key only when
//!   the segment before the version dot equals the requested id, so an
//!   `NM_000088` query still resolves to `NM_000088.3` but `chr1` does
//!   not resolve to `chr1-gene.1`.
//! - exact id matches keep working.
//! - the fallback must not match arbitrary prefixes that stop at a
//!   non-version-boundary character.

#![cfg(feature = "mmap")]

use std::io::Write;

use ferro_hgvs::reference::fasta::MmapFastaProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::ReferenceProvider;

/// Minimal viable FASTA file: one contig, one base. `MmapFastaProvider`
/// needs a real file to memory-map, but the bug under test lives in the
/// transcript registry and never touches the index.
fn tiny_mmap_provider() -> (tempfile::TempDir, MmapFastaProvider) {
    let dir = tempfile::tempdir().expect("create tempdir");
    let fa_path = dir.path().join("tiny.fa");
    let mut f = std::fs::File::create(&fa_path).expect("create fasta");
    writeln!(f, ">chr1").expect("write fasta header");
    writeln!(f, "A").expect("write fasta seq");
    f.sync_all().expect("sync fasta");
    drop(f);
    let provider = MmapFastaProvider::new(&fa_path).expect("mmap fasta");
    (dir, provider)
}

fn synthetic_transcript(id: &str) -> Transcript {
    Transcript::new(
        id.to_string(),
        None,
        Strand::Plus,
        "ACGTACGT".to_string(),
        Some(1),
        Some(8),
        vec![Exon::new(1, 1, 8)],
        Some("chr1".to_string()),
        Some(1),
        Some(8),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    )
}

#[test]
fn get_transcript_does_not_collide_with_chromosome_prefix() {
    let (_dir, mut p) = tiny_mmap_provider();
    p.add_transcript(synthetic_transcript("chr1-gene.1"));
    if let Ok(tx) = p.get_transcript("chr1") {
        panic!(
            "expected ReferenceNotFound for 'chr1'; got transcript {:?}",
            tx.id
        );
    }
}

#[test]
fn get_transcript_resolves_unversioned_query_to_versioned_key() {
    let (_dir, mut p) = tiny_mmap_provider();
    p.add_transcript(synthetic_transcript("NM_000088.3"));
    let tx = p
        .get_transcript("NM_000088")
        .expect("unversioned id should resolve to NM_000088.3");
    assert_eq!(tx.id, "NM_000088.3");
}

#[test]
fn get_transcript_resolves_exact_id() {
    let (_dir, mut p) = tiny_mmap_provider();
    p.add_transcript(synthetic_transcript("NM_000088.3"));
    let tx = p
        .get_transcript("NM_000088.3")
        .expect("exact id should resolve");
    assert_eq!(tx.id, "NM_000088.3");
}

#[test]
fn get_transcript_rejects_non_version_boundary_prefix() {
    let (_dir, mut p) = tiny_mmap_provider();
    p.add_transcript(synthetic_transcript("NM_000088.3"));
    if let Ok(tx) = p.get_transcript("NM_0000") {
        panic!(
            "expected ReferenceNotFound for 'NM_0000'; got transcript {:?}",
            tx.id
        );
    }
}
