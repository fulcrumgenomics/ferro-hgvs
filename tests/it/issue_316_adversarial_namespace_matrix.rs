//! Tests for issue #316 — exercise the new `adversarial_namespace`
//! builder on `tests/common/synthetic.rs`. The builder produces a
//! `MockProvider` whose registered transcript shares its `id` (or the
//! `chromosome` field) with a colliding contig name, with genomic and
//! transcript sequences laid out so a wrong-frame lookup returns
//! observably different bases. This is the failure mode #311 hit and
//! #312 fixed; the matrix here pins post-#312 behavior through the
//! shared builder.
//!
//! Stacks on PR #312. Without #312's `MockProvider::get_sequence` /
//! `get_transcript` fixes, the genomic-coordinate cases here would
//! collapse the merged delins to a single residual SNV.

use crate::common::synthetic::{normalize_to_string, AdversarialNamespaceBuilder};

/// Realistic prefix-collision: a transcript whose id was synthesized
/// from chromosome + gene name (e.g. `chr1-gene.1` from a GFF dump).
#[test]
fn prefix_collision_normalizes_to_canonical_delins() {
    let provider = AdversarialNamespaceBuilder::new("chr1-gene.1", "chr1").build();
    let result = normalize_to_string(provider, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// Already-canonical delins round-trips unchanged under the same
/// prefix-collision fixture.
#[test]
fn prefix_collision_roundtrips_canonical_delins() {
    let provider = AdversarialNamespaceBuilder::new("chr1-gene.1", "chr1").build();
    let result = normalize_to_string(provider, "chr1:g.1000_1002delinsACG");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// Exact-name collision: transcript id is literally the chromosome name.
#[test]
fn exact_name_collision_normalizes_to_canonical_delins() {
    let provider = AdversarialNamespaceBuilder::new("chr1", "chr1").build();
    let result = normalize_to_string(provider, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// Case-sensitivity boundary control: `Chr1-gene.1` does not collide
/// with `chr1` and must keep normalizing correctly.
#[test]
fn case_sensitive_no_collision_still_normalizes() {
    let provider = AdversarialNamespaceBuilder::new("Chr1-gene.1", "chr1").build();
    let result = normalize_to_string(provider, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}

/// No-collision control: a RefSeq-shaped transcript id pinned to the
/// same chromosome must continue to normalize correctly (and is the
/// realistic happy path).
#[test]
fn no_collision_refseq_id_still_normalizes() {
    let provider = AdversarialNamespaceBuilder::new("NM_000001.1", "chr1").build();
    let result = normalize_to_string(provider, "chr1:g.[1000G>A;1001A>C;1002C>G]");
    assert_eq!(result, "chr1:g.1000_1002delinsACG");
}
