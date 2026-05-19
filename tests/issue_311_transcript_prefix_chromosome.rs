//! Tests for issue #311 — `Normalizer.normalize()` drops variants when a
//! transcript `id` is a case-sensitive prefix of the genomic `chromosome`
//! name (e.g. `id = "chr1-gene.1"` with `chromosome = "chr1"`).
//!
//! Root cause: `MockProvider::get_transcript` falls back to a permissive
//! `starts_with(base_id)` lookup intended for unversioned-id resolution
//! (e.g. asking for `NM_000088` and finding `NM_000088.3`). The match has
//! no version-boundary check, so it also matches any transcript whose id
//! literally starts with the requested string — including the chromosome
//! accession when a transcript id happens to share that prefix.
//!
//! When `normalize_genome` calls `provider.get_sequence("chr1", ...)`, the
//! current `MockProvider::get_sequence` tries `get_transcript("chr1")`
//! first. The buggy fallback returns the prefix-colliding transcript, and
//! the call short-circuits to `Transcript::get_sequence(start, end)` —
//! which interprets the genomic coordinates as transcript-relative
//! indices into the transcript's stored sequence. The wrong reference
//! bases drive the post-merge trim down to a single residual SNV.
//!
//! These tests pin the correct end-to-end normalizer behavior and the
//! underlying provider-level lookup contract.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::ReferenceProvider;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Synthetic 1102 bp transcript sequence: 909 bp poly-A leader + 93 bp
/// coding segment + 100 bp poly-T trailer. Positions 1000..=1002 in
/// genomic coordinates (with `genomic_start = 91`) map to transcript
/// positions 910..=912, i.e. `seq[909..912] = "GAC"`. The buggy path
/// reads `seq[999..1002] = "AAG"` instead, which is what produces the
/// erroneous `1001A>C` residual reported in the issue.
const SEGMENT: &str =
    "GACCGGCGGTCTGCAAGTCTCCACCTTCCGAAGCTATCCATAACCGGAACCTATGACCTAAAATCAGTTCTGGGTCAGCTCGGTATCACGAAG";

fn sequence() -> String {
    let mut s = String::with_capacity(1102);
    s.extend(std::iter::repeat_n('A', 909));
    s.push_str(SEGMENT);
    s.extend(std::iter::repeat_n('T', 100));
    assert_eq!(s.len(), 1102);
    assert_eq!(&s[909..912], "GAC");
    s
}

/// Build a `MockProvider` whose only transcript has `id = transcript_id`
/// and `chromosome = chromosome`. Mirrors the issue's reproducer: no
/// top-level `genomic_sequences` entry is registered, so the buggy path
/// in `MockProvider::get_sequence` is the only thing that can resolve a
/// genomic lookup against this provider.
fn provider(transcript_id: &str, chromosome: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    let exons = vec![Exon::with_genomic(1, 1, 1102, 91, 1192)];
    let transcript = Transcript::new(
        transcript_id.to_string(),
        None,
        Strand::Plus,
        sequence(),
        Some(1),
        Some(1102),
        exons,
        Some(chromosome.to_string()),
        Some(91),
        Some(1192),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

fn normalize_with(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

// ===========================================================================
// End-to-end normalizer behavior
// ===========================================================================

/// Transcript id starts with the chromosome name; the normalizer must
/// merge the three adjacent SNVs into a `delins` without dropping
/// positions. Pre-fix: returns `chr1:g.1001A>C`.
#[test]
fn allele_normalizes_when_tx_id_starts_with_chromosome() {
    let p = provider("chr1-gene.1", "chr1");
    let result = normalize_with(p, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// Already-canonical delins on the same genomic span must round-trip
/// unchanged under the same transcript-id collision. Pre-fix: returns
/// `chr1:g.1001A>C`.
#[test]
fn delins_roundtrips_when_tx_id_starts_with_chromosome() {
    let p = provider("chr1-gene.1", "chr1");
    let result = normalize_with(p, "chr1:g.1000_1002delinsACG");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// Negative control: when the transcript id has a typical RefSeq form
/// (no prefix collision with the chromosome name) the same input
/// normalizes correctly today. This pins that the fix does not change
/// behavior for the existing happy path.
#[test]
fn allele_normalizes_when_tx_id_does_not_collide() {
    let p = provider("NM_000001.1", "chr1");
    let result = normalize_with(p, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// Negative control: an exact `chr1` transcript id (no version suffix)
/// is the worst-case collision. Pre-fix: returns `chr1:g.1001A>C`.
#[test]
fn allele_normalizes_when_tx_id_equals_chromosome() {
    let p = provider("chr1", "chr1");
    let result = normalize_with(p, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// Sanity check on the case-sensitivity boundary called out in the
/// issue: `Chr1-gene.1` (capital C) does not collide with `chr1` today
/// and must continue not to collide after the fix.
#[test]
fn allele_normalizes_when_tx_id_differs_only_in_case() {
    let p = provider("Chr1-gene.1", "chr1");
    let result = normalize_with(p, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

// ===========================================================================
// Provider-level lookup contract
// ===========================================================================

/// Pin the exact provider-level invariant: `get_transcript("chr1")` must
/// not resolve to `chr1-gene.1` via the unversioned-prefix fallback.
/// The fallback is meant to bridge `NM_000088` ↔ `NM_000088.3`, where
/// the only character following the requested base id is a version dot.
#[test]
fn get_transcript_does_not_collide_with_chromosome_prefix() {
    let p = provider("chr1-gene.1", "chr1");
    if let Ok(tx) = p.get_transcript("chr1") {
        panic!(
            "expected ReferenceNotFound for 'chr1'; got transcript {:?}",
            tx.id
        );
    }
}

/// The unversioned-prefix fallback must still resolve a versionless
/// query to a versioned key. This is the case the fallback was added
/// for in the first place and must continue to work.
#[test]
fn get_transcript_resolves_unversioned_query_to_versioned_key() {
    let p = provider("NM_000088.3", "chr1");
    let tx = p
        .get_transcript("NM_000088")
        .expect("unversioned id should resolve to NM_000088.3");
    assert_eq!(tx.id, "NM_000088.3");
}

/// An exact id match must still resolve. This guards against an
/// over-strict fix that breaks the primary lookup path.
#[test]
fn get_transcript_resolves_exact_id() {
    let p = provider("NM_000088.3", "chr1");
    let tx = p
        .get_transcript("NM_000088.3")
        .expect("exact id should resolve");
    assert_eq!(tx.id, "NM_000088.3");
}

/// A query that is a literal prefix of a versioned id but stops at a
/// non-version-boundary character (e.g. asking for `NM_0000` against
/// `NM_000088.3`) must NOT resolve. The unversioned-prefix fallback is
/// strictly about the version-suffix boundary, not arbitrary prefixes.
#[test]
fn get_transcript_rejects_non_version_boundary_prefix() {
    let p = provider("NM_000088.3", "chr1");
    if let Ok(tx) = p.get_transcript("NM_0000") {
        panic!(
            "expected ReferenceNotFound for 'NM_0000'; got transcript {:?}",
            tx.id
        );
    }
}

/// A versioned query that does not exactly match any stored key must
/// NOT cross-version fall back. The base-id fallback exists only to
/// resolve unversioned queries (e.g. `NM_000088`) against a versioned
/// store; a versioned query (`NM_000088.5`) that misses an exact match
/// must error rather than silently return a different version
/// (`NM_000088.3`).
#[test]
fn get_transcript_rejects_versioned_query_against_other_version() {
    let p = provider("NM_000088.3", "chr1");
    if let Ok(tx) = p.get_transcript("NM_000088.5") {
        panic!(
            "expected ReferenceNotFound for 'NM_000088.5'; got transcript {:?}",
            tx.id
        );
    }
}

/// Same contract for `get_protein_sequence`: an unversioned query must
/// still resolve to a stored versioned key, but a versioned query that
/// misses an exact match must NOT cross-version fall back.
#[test]
fn get_protein_sequence_versioning_contract() {
    let mut p = MockProvider::new();
    p.add_protein("NP_000079.2", "MFSFVDLRLLLLLAATALLTHGQEEGQVEGQDEDIPP");

    let seq = p
        .get_protein_sequence("NP_000079", 0, 5)
        .expect("unversioned accession should resolve to NP_000079.2");
    assert_eq!(seq, "MFSFV");

    let seq = p
        .get_protein_sequence("NP_000079.2", 0, 5)
        .expect("exact accession should resolve");
    assert_eq!(seq, "MFSFV");

    if let Ok(seq) = p.get_protein_sequence("NP_000079.5", 0, 5) {
        panic!(
            "expected ProteinReferenceNotAvailable for 'NP_000079.5'; got {:?}",
            seq
        );
    }
}
