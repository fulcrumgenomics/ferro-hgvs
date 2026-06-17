//! Tests for issue #315 — `FastaProvider::get_sequence` routes every
//! lookup through the transcript registry before consulting the FASTA
//! index. Exact-match only (no `starts_with` fallback), but if any
//! caller registers a transcript with `id = "chrN"` alongside a FASTA
//! that also contains `chrN`, the transcript wins and a genomic
//! coordinate lookup silently reads transcript-relative bases — the
//! same wrong-coordinate-frame bug #311 hit on `MockProvider`.
//!
//! These tests pin the corrected dispatch contract:
//!
//! - when the id is a known contig (in the FASTA `.fai` index, or as the
//!   `chromosome` field on any registered transcript), `get_sequence`
//!   resolves strictly against the FASTA path.
//! - transcript-id lookups for ids that are NOT known contigs keep
//!   returning the transcript's cached sequence (no regression).
//! - a known contig with no FASTA backing surfaces as `Err`, which
//!   `normalize_genome` already falls back to canonicalize-only.

use std::io::Write;

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::{FastaProvider, ReferenceProvider};

/// FASTA file with a single contig `chr1` whose sequence is 12 bytes
/// of `A`. `MmapFastaProvider` and `FastaProvider` both read from a
/// real file; `tempfile::TempDir` keeps the file alive for the
/// duration of the test.
fn provider_with_chr1_fasta() -> (tempfile::TempDir, FastaProvider) {
    let dir = tempfile::tempdir().expect("create tempdir");
    let fa_path = dir.path().join("tiny.fa");
    let mut f = std::fs::File::create(&fa_path).expect("create fasta");
    writeln!(f, ">chr1").expect("write fasta header");
    writeln!(f, "AAAAAAAAAAAA").expect("write fasta seq");
    f.sync_all().expect("sync fasta");
    drop(f);
    let provider = FastaProvider::new(&fa_path).expect("load fasta");
    (dir, provider)
}

fn provider_no_chr1_fasta() -> (tempfile::TempDir, FastaProvider) {
    let dir = tempfile::tempdir().expect("create tempdir");
    let fa_path = dir.path().join("tiny.fa");
    let mut f = std::fs::File::create(&fa_path).expect("create fasta");
    writeln!(f, ">chr2").expect("write fasta header");
    writeln!(f, "GGGGGGGG").expect("write fasta seq");
    f.sync_all().expect("sync fasta");
    drop(f);
    let provider = FastaProvider::new(&fa_path).expect("load fasta");
    (dir, provider)
}

fn synthetic_transcript(id: &str, chromosome: Option<&str>) -> Transcript {
    Transcript::new(
        id.to_string(),
        None,
        Strand::Plus,
        "TTTTTTTT".to_string(),
        Some(1),
        Some(8),
        vec![Exon::new(1, 1, 8)],
        chromosome.map(str::to_string),
        Some(1),
        Some(8),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    )
}

/// Exact id collision: a transcript with id `chr1` is registered on
/// the same provider whose FASTA index also has `chr1`. The FASTA
/// bytes (`AAA`) must win over the transcript's cached bases (`TTT`).
#[test]
fn get_sequence_prefers_fasta_index_over_colliding_transcript_id() {
    let (_dir, mut p) = provider_with_chr1_fasta();
    p.add_transcript(synthetic_transcript("chr1", Some("chr1")));
    let seq = p
        .get_sequence("chr1", 0, 3)
        .expect("chr1 lookup should resolve via FASTA");
    assert_eq!(seq, "AAA", "expected FASTA bytes, got transcript bases");
}

/// Soft-collision: a transcript declares `chr1` as its `chromosome`
/// field, but the FASTA index has no `chr1` entry. A `get_sequence`
/// call on `chr1` must NOT silently slide into the transcript's bases
/// for an unrelated transcript id — it must Err so callers can fall
/// back to canonicalize-only.
#[test]
fn get_sequence_for_known_contig_without_fasta_backing_errors() {
    let (_dir, mut p) = provider_no_chr1_fasta();
    p.add_transcript(synthetic_transcript("NM_000088.3", Some("chr1")));
    if let Ok(seq) = p.get_sequence("chr1", 0, 3) {
        panic!("expected ReferenceNotFound for 'chr1'; got {:?}", seq);
    }
}

/// Negative control: a transcript id with no colliding chromosome must
/// keep returning the transcript's cached sequence.
#[test]
fn get_sequence_for_transcript_id_returns_transcript_bases() {
    let (_dir, mut p) = provider_with_chr1_fasta();
    p.add_transcript(synthetic_transcript("NM_000088.3", Some("chr1")));
    let seq = p
        .get_sequence("NM_000088.3", 0, 3)
        .expect("NM_000088.3 lookup should resolve via transcript registry");
    assert_eq!(seq, "TTT");
}

/// Negative control: a FASTA-index lookup with no transcript
/// registered for that name keeps returning FASTA bytes.
#[test]
fn get_sequence_for_fasta_only_contig_returns_fasta_bases() {
    let (_dir, p) = provider_with_chr1_fasta();
    let seq = p
        .get_sequence("chr1", 0, 3)
        .expect("chr1 lookup should resolve via FASTA");
    assert_eq!(seq, "AAA");
}
