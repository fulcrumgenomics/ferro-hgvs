//! #2155 — the derivation surface (`from_sequences`) must collapse a
//! payload-coincidence change the same way the normalization surface
//! (`Normalizer::normalize`) does, on every DNA axis it can reach (`g.`/`m.`).
//!
//! `from_sequences` never runs `canonicalize_from_sequence`'s coalesce
//! passes — it reaches `derive_block_members` instead, which used to stop at
//! `shift_pieces`/`coalesce_adjacent_pieces`/`shrink_pieces_to_differences` and
//! emit the fragmented split. This is the same block
//! `issue_2155_payload_coincidence_all_dna.rs` proves collapses on `normalize`;
//! this file proves it now also collapses on `from_sequences`, and that the two
//! surfaces agree on the string.
//!
//! Block at 1-based 10_17: `CTTAGTTA -> AAACAAAC` (equal length, a payload
//! coincidence per `DNA/delins.md:44-47`).

use ferro_hgvs::{from_sequences, parse_hgvs, FromSequencesOptions, MockProvider, Normalizer};

const REF: &str = "AGAACCCCCCTTAGTTAAGAACAAAAGCAACAATCTTCGTGGTCCTGG";
const ALT: &str = "AGAACCCCCAAACAAACAGAACAAAAGCAACAATCTTCGTGGTCCTGG";

/// The genomic axis: `from_sequences` on an arbitrary genomic accession must
/// collapse to the single spanning `delins`, not the fragmented split
/// (`g.[10_12delinsAA;14_16delinsCAA;17_18insC]`).
#[test]
fn the_genomic_axis_collapses_to_the_single_spanning_delins() {
    let derived = from_sequences("TEMPLATE", 1, REF, ALT, &FromSequencesOptions::default())
        .expect("from_sequences must derive this block");
    assert_eq!(derived.to_string(), "TEMPLATE:g.10_17delinsAAACAAAC");
}

/// The mitochondrial axis, on `from_sequences`'s own `m.` accession
/// (`Accession::is_mitochondrial` -- `NC_012920`/`NC_001807` only).
#[test]
fn the_mitochondrial_axis_collapses_to_the_single_spanning_delins() {
    let derived = from_sequences("NC_012920.1", 1, REF, ALT, &FromSequencesOptions::default())
        .expect("from_sequences must derive this block");
    assert_eq!(derived.to_string(), "NC_012920.1:m.10_17delinsAAACAAAC");
}

/// Convergence: `from_sequences`'s derivation of the block and
/// `Normalizer::normalize`'s re-derivation of the single spanning `delins` over
/// the same reference must agree, byte for byte. If they diverge, the coalesce
/// order/subset `derive_block_members` applies is not the one
/// `canonicalize_from_sequence` applies.
#[test]
fn from_sequences_converges_with_normalize_on_the_same_block() {
    let derived = from_sequences("TEMPLATE", 1, REF, ALT, &FromSequencesOptions::default())
        .expect("from_sequences must derive this block")
        .to_string();

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", REF.to_string());
    let input = parse_hgvs("TEMPLATE:g.10_17delinsAAACAAAC").expect("parses");
    let normalizer = Normalizer::new(provider);
    let normalized = normalizer
        .normalize(&input)
        .expect("normalizes")
        .to_string();

    // Feeding the collapsed form only pins idempotency of the spanning delins.
    // The fragmented spelling this pass exists to collapse must reach the same
    // string on `normalize`, or a normalize-side regression would slip past this
    // convergence check while `issue_2155_payload_coincidence_all_dna.rs` carried
    // the whole burden.
    let fragmented =
        parse_hgvs("TEMPLATE:g.[10_12delinsAA;14_16delinsCAA;17_18insC]").expect("parses");
    assert_eq!(
        normalizer
            .normalize(&fragmented)
            .expect("normalizes")
            .to_string(),
        "TEMPLATE:g.10_17delinsAAACAAAC",
        "normalize must collapse the fragmented spelling to the spanning delins",
    );

    assert_eq!(
        derived, normalized,
        "from_sequences and normalize must converge on the same block"
    );
    assert_eq!(derived, "TEMPLATE:g.10_17delinsAAACAAAC");
}
